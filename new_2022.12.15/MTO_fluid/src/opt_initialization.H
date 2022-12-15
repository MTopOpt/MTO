Info << "Initialization\n"<< endl;
std::remove("Voluse.txt");
std::remove("Time.txt");
std::remove("MeanT.txt");
std::remove("maxT.txt");
std::remove("PowerDiss.txt");

int i = 0, opt = 0;
double del = 0;
scalar n(0), N(0), eta0, eta5, eta1, y0, y5;
n = x.size(); // n:length of x in the i-th core
N = x.size(); // N:design freedom
reduce(N, sumOp<scalar>());

scalar C(0), V(0), area = gSum(mesh.V()), MeanT(0), PowerDiss(0), maxT(0), sumT(0);
double len = std::pow(area / N, 1.0 / 3);
b.value() = 1.0 / (filterR * len / 3.464) / (filterR * len / 3.464);
//***********************************//MMA
int m = 1; // number of constraint
MMA mma(n, m);
mma.asyminit = mma_init;
mma.asymdec = mma_dec;
mma.asyminc = mma_inc;
mma.raa0 = raa0;
mma.movelimit = movlim;
std::vector<double> xmma(n);
std::vector<double> dfdx(n);
std::vector<double> gx(m);
std::vector<std::vector<double>> dgdx(m);
for (int i = 0; i < m; i++)
{
  dgdx[i].resize(n);
}
