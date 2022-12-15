#include "MMA.h"

void MMA::MMAsolver(std::vector<double> &xval,
              std::vector<double> &dfdx,
              std::vector<double> &g,
              std::vector<std::vector<double>> &dgdx)
{
    if (iter==1)
    {
        xold1=xval;
        xold2=xval;
    }
    
    Update(xval, dfdx, g, dgdx);
    xold2=xold1;
    xold1=xval;
    iter++;
}

MMA::MMA(int NvarLocal, int Mcons)
    : n(NvarLocal), m(Mcons), iter(1), y(m), lam(m), mu(m),
      s(2 * m), low(n), upp(n), xmax(n, 1), xmin(n, 0), alpha(n), beta(n),
      p0(n), q0(n), pij(n * m), qij(n * m), b(m),
      grad(m), hess(m * m), a(m), c(m),xold1(n),xold2(n)
{
    int NvarGlab = 0;
    MPI_Allreduce(&n, &NvarGlab, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    epsimin = sqrt(NvarGlab + m) * 1e-9;
    for (int i = 0; i < Mcons; i++)
    {
        // a0[i] = 1.0;
        a[i] = 0;
        c[i] = 1000.0;
        // d[i]  = 1;
    }
}

void MMA::Update(std::vector<double> &xval,
              std::vector<double> &dfdx,
              std::vector<double> &g,
              std::vector<std::vector<double>> &dgdx)
{
    // Generate the subproblem
    GenSub(xval, dfdx, g, dgdx);

    // Solve the dual with an interior point method
    SolveDIP(xval);

}

void MMA::SolveDIP(std::vector<double> &x)
{
    for (int j = 0; j < m; j++)
    {
        lam[j] = c[j] / 2.0;
        mu[j] = 1.0;
    }

    const double tol = epsimin; // 1.0e-9*sqrt(m+n);
    double epsi = 1.0;
    double err = 1.0;
    int loop;

    while (epsi > tol)
    {
        loop = 0;
        while (err > 0.9 * epsi && loop < 100)
        {
            loop++;

            // Set up Newton system
            XYZofLAMBDA(x);
            DualGrad(x);
            for (int j = 0; j < m; j++)
            {
                grad[j] = -1.0 * grad[j] - epsi / lam[j];
            }
            DualHess(x);

            // Solve Newton system
            if (m > 1.5)
            {
                Factorize(hess.data(), m);
                Solve(hess.data(), grad.data(), m);
                for (int j = 0; j < m; j++)
                {
                    s[j] = grad[j];
                }
            }
            else if (m > 0)
            {
                s[0] = grad[0] / hess[0];
            }

            // Get the full search direction
            for (int i = 0; i < m; i++)
            {
                s[m + i] = -mu[i] + epsi / lam[i] - s[i] * mu[i] / lam[i];
            }

            // Perform linesearch and update lam and mu
            DualLineSearch();

            XYZofLAMBDA(x);

            // Compute KKT res
            err = DualResidual(x, epsi);
        }
        epsi = epsi * 0.1;
    }
}
double MMA::DualResidual(std::vector<double> &x, double epsi)
{

    double *res = new double[2 * m];
    double *resGlb = new double[2 * m];

    for (int j = 0; j < m; j++)
    {
        res[j] = 0;
        res[j + m] = 0;
        for (int i = 0; i < n; i++)
        {
            res[j] += pij[i * m + j] / (upp[i] - x[i]) + qij[i * m + j] / (x[i] - low[i]);
        }
    }

    MPI_Allreduce(res, resGlb, 2 * m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for (int i = 0; i < 2 * m; i++)
    {
        res[i] = resGlb[i];
    }
	for (int j = 0; j < m; j++)
    {
        res[j] += -b[j] - a[j] * z - y[j] + mu[j];
        res[j + m] += mu[j] * lam[j] - epsi;
    }
    double nrI = 0.0;
    for (int i = 0; i < 2 * m; i++)
    {
        if (nrI < std::abs(res[i]))
        {
            nrI = std::abs(res[i]);
        }
    }

    delete[] res;
    delete[] resGlb;

    return nrI;
}
void MMA::DualLineSearch()
{

    double theta = 1.005;
    for (int i = 0; i < m; i++)
    {
        if (theta < -1.01 * s[i] / lam[i])
        {
            theta = -1.01 * s[i] / lam[i];
        }
        if (theta < -1.01 * s[i + m] / mu[i])
        {
            theta = -1.01 * s[i + m] / mu[i];
        }
    }
    theta = 1.0 / theta;

    for (int i = 0; i < m; i++)
    {
        lam[i] = lam[i] + theta * s[i];
        mu[i] = mu[i] + theta * s[i + m];
    }
}
void MMA::DualHess(std::vector<double> &x)
{

    double *df2 = new double[n];
    double *PQ = new double[n * m];

    for (int i = 0; i < n; i++)
    {
        double pjlam = p0[i];
        double qjlam = q0[i];
        for (int j = 0; j < m; j++)
        {
            pjlam += pij[i * m + j] * lam[j];
            qjlam += qij[i * m + j] * lam[j];
            PQ[i * m + j] = pij[i * m + j] / pow(upp[i] - x[i], 2.0) - qij[i * m + j] / pow(x[i] - low[i], 2.0);
        }
        df2[i] = -1.0 / (2.0 * pjlam / pow(upp[i] - x[i], 3.0) + 2.0 * qjlam / pow(x[i] - low[i], 3.0));
        double xp = (sqrt(pjlam) * low[i] + sqrt(qjlam) * upp[i]) / (sqrt(pjlam) + sqrt(qjlam));
        if (xp < alpha[i])
        {
            df2[i] = 0.0;
        }
        if (xp > beta[i])
        {
            df2[i] = 0.0;
        }
    }

    // Create the matrix/matrix/matrix product: PQ^T * diag(df2) * PQ
    double *tmp = new double[n * m];
    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {
            tmp[j * n + i] = 0.0;
            tmp[j * n + i] += PQ[i * m + j] * df2[i];
        }
    }

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            hess[i * m + j] = 0.0;
            for (int k = 0; k < n; k++)
            {
                hess[i * m + j] += tmp[i * n + k] * PQ[k * m + j];
            }
        }
    }

    double lamai = 0.0;
    std::vector<double> hessGlb(hess);
    MPI_Allreduce(&hess.front(), &hessGlb.front(), m * m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for (int i = 0; i < m * m; i++)
    {
        hess[i] = hessGlb[i];
    }
    for (int j = 0; j < m; j++)
    {
        if (lam[j] < 0.0)
        {
            lam[j] = 0.0;
        }
        lamai += lam[j] * a[j];
        if (lam[j] > c[j])
        {
            hess[j * m + j] += -1.0;
        }
        hess[j * m + j] += -mu[j] / lam[j];
    }

    if (lamai > 0.0)
    {
        for (int j = 0; j < m; j++)
        {
            for (int k = 0; k < m; k++)
            {
                hess[j * m + k] += -10.0 * a[j] * a[k];
            }
        }
    }

    // pos def check
    double HessTrace = 0.0;
    for (int i = 0; i < m; i++)
    {
        HessTrace += hess[i * m + i];
    }
    double HessCorr = 1e-4 * HessTrace / m;

    if (-1.0 * HessCorr < 1.0e-7)
    {
        HessCorr = -1.0e-7;
    }

    for (int i = 0; i < m; i++)
    {
        hess[i * m + i] += HessCorr;
    }

    delete[] df2;
    delete[] PQ;
    delete[] tmp;
}
void MMA::Factorize(double *K, int n)
{

    for (int s = 0; s < n - 1; s++)
    {
        for (int i = s + 1; i < n; i++)
        {
            K[i * n + s] = K[i * n + s] / K[s * n + s];
            for (int j = s + 1; j < n; j++)
            {
                K[i * n + j] = K[i * n + j] - K[i * n + s] * K[s * n + j];
            }
        }
    }
}

void MMA::Solve(double *K, double *x, int m)
{

    for (int i = 1; i < m; i++)
    {
        double a = 0.0;
        for (int j = 0; j < i; j++)
        {
            a = a - K[i * m + j] * x[j];
        }
        x[i] = x[i] + a;
    }

    x[m - 1] = x[m- 1] / K[(m - 1) * m + (m - 1)];
    for (int i = m - 2; i >= 0; i--)
    {
        double a = x[i];
        for (int j = i + 1; j < m; j++)
        {
            a = a - K[i * m + j] * x[j];
        }
        x[i] = a / K[i * m + i];
    }
}



void MMA::DualGrad(std::vector<double> &x)
{
    for (int j = 0; j < m; j++)
    {
        grad[j] =0;
        for (int i = 0; i < n; i++)
        {
            grad[j] += pij[i * m + j] / (upp[i] - x[i]) + qij[i * m + j] / (x[i] - low[i]);
        }
    }
    std::vector<double> gradGlb(grad);
    MPI_Allreduce(&grad.front(), &gradGlb.front(), m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int j = 0; j < m; j++)
    {
        grad[j] = gradGlb[j] - b[j] - a[j] * z - y[j];
    }
}

void MMA::XYZofLAMBDA(std::vector<double> &x)
{
    double lamai = 0.0;
    for (int i = 0; i < m; i++)
    {
        if (lam[i] < 0.0)
        {
            lam[i] = 0;
        }
        y[i] = std::max(0.0, lam[i] - c[i]); // Note y=(lam-c)/d - however d is fixed at one !!
        lamai += lam[i] * a[i];
    }
    z = std::max(0.0, 10.0 * (lamai - 1.0)); // SINCE a0 = 1.0

    for (int i = 0; i < n; i++)
    {
        double pjlam = p0[i];
        double qjlam = q0[i];
        for (int j = 0; j < m; j++)
        {
            pjlam += pij[i * m + j] * lam[j];
            qjlam += qij[i * m + j] * lam[j];
        }
        x[i] = (sqrt(pjlam) * low[i] + sqrt(qjlam) * upp[i]) / (sqrt(pjlam) + sqrt(qjlam));
        if (x[i] < alpha[i])
        {
            x[i] = alpha[i];
        }
        if (x[i] > beta[i])
        {
            x[i] = beta[i];
        }
    }
}

void MMA::GenSub(std::vector<double> &xval, std::vector<double> &dfdx,
                  std::vector<double> &g, std::vector<std::vector<double>> &dgdx)
{
    for (int i = 0; i < n; i++)
    {
        xmax[i] = std::min(1., xval[i] + movelimit);
        xmin[i] = std::max(0., xval[i] - movelimit);
    }
    // Set asymptotes
    if (iter < 3)
    {
        for (int i = 0; i < n; i++)
        {
            low[i] = xval[i] - asyminit * (xmax[i] - xmin[i]);
            upp[i] = xval[i] + asyminit * (xmax[i] - xmin[i]);
        }
    }
    else
    {
        for (int i = 0; i < n; i++)
        {
            double zzz = (xval[i] - xold1[i]) * (xold1[i] - xold2[i]);
            double gamma;
            if (zzz < 0.0)
            {
                gamma = asymdec;
            }
            else if (zzz > 0.0)
            {
                gamma = asyminc;
            }
            else
            {
                gamma = 1.0;
            }
            low[i] = xval[i] - gamma * (xold1[i] - low[i]);
            upp[i] = xval[i] + gamma * (upp[i] - xold1[i]);

            double xmami = std::max(xmamieps, xmax[i] - xmin[i]);
            if (RobustAsymptotesType < 0.5)
            {
                low[i] = std::max(low[i], xval[i] - 10.0 * xmami);
                low[i] = std::min(low[i], xval[i] - 0.01 * xmami);
                upp[i] = std::min(upp[i], xval[i] + 10.0 * xmami);
                upp[i] = std::max(upp[i], xval[i] + 0.01 * xmami);
            }
            else // if (RobustAsymptotesType == 1)
            {
                low[i] = std::max(low[i], xval[i] - 100.0 * xmami);
                low[i] = std::min(low[i], xval[i] - 1.0e-4 * xmami);
                upp[i] = std::max(upp[i], xval[i] + 1.0e-4 * xmami);
                upp[i] = std::min(upp[i], xval[i] + 100.0 * xmami);
                double xmi = xmin[i] - 1.0e-5;
                double xma = xmax[i] + 1.0e-5;
                if (xval[i] < xmi)
                {
                    low[i] = xval[i] - (xma - xval[i]) / 0.9;
                    upp[i] = xval[i] + (xma - xval[i]) / 0.9;
                }
                if (xval[i] > xma)
                {
                    low[i] = xval[i] - (xval[i] - xmi) / 0.9;
                    upp[i] = xval[i] + (xval[i] - xmi) / 0.9;
                }
            }
        }
    }

    // Set bounds and the coefficients for the approximation
    // double raa0 = 0.5*1e-6;
    for (int i = 0; i < n; i++)
    {
        // Compute bounds alpha and beta
        alpha[i] = std::max(xmin[i], low[i] + albefa * (xval[i] - low[i]));
        // alpha[i] = std::max(alpha[i], xval[i] - movelimit * (xmax[i] - xmin[i]));
        // alpha[i] = std::min(alpha[i], xmax[i]);

        beta[i] = std::min(xmax[i], upp[i] - albefa * (upp[i] - xval[i]));
        // beta[i] = std::min(beta[i], xval[i] + movelimit * (xmax[i] - xmin[i]));
        // beta[i] = std::max(beta[i], xmin[i]);

        // Objective function
        {
            double dfdxp = std::max(0.0, dfdx[i]);
            double dfdxm = std::max(0.0, -1.0 * dfdx[i]);
            double xmamiinv = 1.0 / std::max(xmamieps, xmax[i] - xmin[i]);
            double pq = 0.001 * std::abs(dfdx[i]) + raa0 * xmamiinv;
            // double pq = 0.001 * std::abs(dfdx[i]) + raa0 ;
            p0[i] = std::pow(upp[i] - xval[i], 2.0) * (dfdxp + pq);
            q0[i] = std::pow(xval[i] - low[i], 2.0) * (dfdxm + pq);
        }

        // Constraints
        for (int j = 0; j < m; j++)
        {
            double dgdxp = std::max(0.0, dgdx[j][i]);
            double dgdxm = std::max(0.0, -1.0 * dgdx[j][i]);
            double xmamiinv = 1.0 / std::max(xmamieps, xmax[i] - xmin[i]);
            double pq = 0.001 * std::abs(dgdx[j][i]) + raa0 * xmamiinv;
            // double pq = 0.001 * std::abs(dgdx[i * m + j]) + raa0 ;
            pij[i * m + j] = std::pow(upp[i] - xval[i], 2.0) * (dgdxp + pq);
            qij[i * m + j] = std::pow(xval[i] - low[i], 2.0) * (dgdxm + pq);
        }
    }

    // The constant for the constraints
    for (int j = 0; j < m; j++)
    {
        b[j] = 0;
        for (int i = 0; i < n; i++)
        {
            b[j] += pij[i * m + j] / (upp[i] - xval[i]) + qij[i * m + j] / (xval[i] - low[i]);
        }
    }

    std::vector<double> sumb(b);
    MPI_Allreduce(&b.front(), &sumb.front(), m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for (int j = 0; j < m; j++)
    {
        b[j] = sumb[j] - g[j];
    }
}