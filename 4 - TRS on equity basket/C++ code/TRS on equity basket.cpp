#include <algorithm>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <math.h>
#include <any>
#include <Eigen/Cholesky>
#include <EigenRand/EigenRand>
#include <numeric>

typedef std::vector<std::unordered_map<std::string, std::any>> typePaths;

//Market datas
auto M = 3;

//TRS datas
auto T = 1.0;
std::vector<double> TRSpaymentSchedule = { 0,25, 50, 75, 99 };
auto sTRS = 0.01;

//Parameters
auto deltaT = 0.01;
auto lT = static_cast<int>(T / deltaT);
auto nbSimul = 100;

typePaths pathsGeneration(int nbSimul, int lT, double r0, const std::vector<double>& kappa,
    const std::vector<double>& theta, const std::vector<double>& lamb, double nu, double etha,
    double alpha, const Eigen::MatrixXd& L, const std::vector<double>& S0, const std::vector<double>& sigma0,
    std::vector<double>& rhoSsigma,const std::vector<double>& rhoSr)
{
    typePaths paths;
    for (auto i = 0; i != nbSimul; i++)
    {
        std::default_random_engine generator;
        Eigen::MatrixXd x = Eigen::Rand::normal<Eigen::MatrixXd>(M, lT, generator, 0.0, 1.0);
        Eigen::MatrixXd phiS = L.cast<double>() * x.cast<double>();
        Eigen::MatrixXd phitild = Eigen::Rand::normal<Eigen::MatrixXd>(M, lT, generator, 0.0, 1.0);
        std::vector<double> phihat;
        for (auto k = 0; k != lT; k++)
        {
            phihat.emplace_back(std::normal_distribution<double>{0, 1}(generator));
        }
        Eigen::MatrixXf S(M,lT);
        Eigen::MatrixXf sigma(M,lT);
        std::vector<double> r(lT);
        std::vector<double> intr(lT);
        for (auto k = 0; k != M;k++)
        {
            S(k, 0) = S0[k];
            sigma(k, 0) = sigma0[k];
        }
        r[0] = r0;
        intr[0] = r[0];
        for (auto k = 0; k != (lT - 1); k++)
        {
            for (auto j = 0; j != M; j++)
            {
                r[k + 1] = (nu - etha * r[k]) * deltaT + sqrt(alpha * r[k]) * (rhoSr[j] * sqrt(deltaT) * phiS(j,k) + 
                    sqrt(1 - rhoSr[j] * rhoSr[j]) * sqrt(deltaT) * phihat[k]) + r[k];
                intr[k + 1] = intr[k] + r[k + 1];
                sigma(j,k+1) = kappa[j] * (theta[j] - sigma(j,k)) * deltaT + lamb[j] * sqrt(sigma(j,k)) * (rhoSsigma[j] * sqrt(deltaT) * phiS(j,k) +
                    sqrt(1 - rhoSsigma[j] * rhoSsigma[j]) * sqrt(deltaT) * phitild(j,k)) + sigma(j,k);
                S(j,k+1) = S(j,k) * (r[k] * deltaT + sqrt(sigma(j,k) * deltaT) * phiS(j,k)) + S(j,k);
            }
        }
        std::unordered_map<std::string, std::any> interRes = { {"Asset",S},{"Volatility",sigma},{"Interest rate",r},{"Cum. interest rate",intr} };
        paths.emplace_back(interRes);
    }
    return paths;
}

std::vector<double> alphaBetaTRS(int TRSpaydate, int prevTRSpaydate, const typePaths& paths, int nbSimul)
{
    auto alphak = 0.0;
    auto betak = 0.0;
    for (auto i = 0; i != nbSimul; i++)
    {
        auto sum1 = 0.0;
        auto sum2 = 0.0;
        for (auto j = 0; j != M; j++)
        {
            sum1 += std::any_cast<Eigen::MatrixXf>(paths[i].at("Asset"))(j,TRSpaydate);
            sum2 += std::any_cast<Eigen::MatrixXf>(paths[i].at("Asset"))(j, prevTRSpaydate);
        }
        auto res1 = (sum1 - sum2) / (exp(std::any_cast<std::vector<double>>(paths[i].at("Cum. interest rate"))[TRSpaydate]));
        auto res2 = (sum2) / (exp(std::any_cast<std::vector<double>>(paths[i].at("Cum. interest rate"))[TRSpaydate]));
        alphak += res1;
        betak += res2;
    }
    alphak /= nbSimul;
    betak /= nbSimul;
    std::vector<double> res;
    res.emplace_back(alphak);
    res.emplace_back(betak);
    return res;
}

int getIndex(const std::vector<double>& v, int K)
{
    auto it = find(v.begin(), v.end(), K);
    if (it != v.end())
    {
        int index = it - v.begin();
        return index;
    }
    else
    {
        return 1;
    }
}

std::vector<double> TRS(const std::vector<double>& TRSpaymentSchedule, double sTRS, const typePaths& paths, int nbSimul)
{
    auto alpha = 0.0;
    auto beta = 0.0;
    for (auto s = TRSpaymentSchedule.begin(); s != TRSpaymentSchedule.end() - 1; s++)
    {
        auto payDate = getIndex(TRSpaymentSchedule, *s) + 1;
        auto prevPayDate = getIndex(TRSpaymentSchedule, *s);
        auto abRes = alphaBetaTRS(TRSpaymentSchedule[payDate], TRSpaymentSchedule[prevPayDate], paths, nbSimul);
        alpha += abRes[0];
        beta += abRes[1] * (TRSpaymentSchedule[payDate] - TRSpaymentSchedule[prevPayDate]);
    }
    auto resTemp = alpha - sTRS * beta;
    auto impliedRate = alpha / beta;
    std::vector<double> res;
    res.emplace_back(resTemp);
    res.emplace_back(impliedRate);
    return res;
}

int main()
{
    //Market datas
    std::vector<double> S0(M);
    std::iota(S0.begin(), S0.end(), 100);
    std::vector<double> sigma0(M);
    std::fill(sigma0.begin(), sigma0.end(), 0.1);
    auto r0 = 0.05;

    //Model datas
    std::vector<double> kappa(M);
    std::vector<double> theta(M);
    std::vector<double> lamb(M);
    std::fill(kappa.begin(), kappa.end(), 0.1);
    std::fill(theta.begin(), theta.end(), 0.1);
    std::fill(lamb.begin(), lamb.end(), 0.1);
    auto nu = 0.1;
    auto etha = 0.1;
    auto alpha = 0.02;

    //Model correlations"
    Eigen::Matrix3f corrSS;
    corrSS << 1, 0.2, 0.1,
              0.2, 1, 0.3,
              0.1, 0.3, 1;
    Eigen::MatrixXd L(corrSS.cast<double>().llt().matrixL());
    std::vector<double> rhoSsigma(M);
    std::vector<double> rhoSr(M);
    std::fill(rhoSsigma.begin(), rhoSsigma.end(), 0.1);
    std::fill(rhoSr.begin(), rhoSr.end(), 0.1);

    //Pricing
    auto paths = pathsGeneration(nbSimul, lT, r0, kappa,
        theta, lamb, nu, etha,
        alpha, L, S0, sigma0,
        rhoSsigma, rhoSr);
    auto TRSres = TRS(TRSpaymentSchedule, sTRS, paths, nbSimul);
    std::cout << "TRS market value is" << TRSres[0] << std::endl;
    std::cout << "TRS implied rate is" << TRSres[1] << std::endl;
}