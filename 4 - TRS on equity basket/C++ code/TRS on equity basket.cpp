#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <math.h>

typedef std::vector<std::unordered_map<std::string, std::any>> typePaths;

//Market datas
auto S0 = 100.0;
auto sigma0 = 0.1;
r0 = 0.05;

//TRS datas
auto T = 1.0;
auto TRSpaymentSchedule = { 25, 50, 75, 99 };
auto sTRS = 0.01;

//Model datas
auto kappa = 0.1;
auto theta = 0.1;
auto lamb = 0.1;
auto nu = 0.1;
auto etha = 0.1;
auto alpha = 0.02;
auto rhoSsigma = 0.1;
auto rhoSr = 0.1;
auto deltaT = 0.01;
auto lT = int(T / deltaT);
auto nbSimul = 1000;

typePaths pathsGeneration(int nbSimul, int lT, double S0, double sigma0, double r0, double kappa, double theta, double lamb, double nu, double etha, double alpha, double rhoSsigma, double rhoSr)
{
    typePaths paths;
    for (i = 0; i != nbSimul; i++)
    {
        std::vector<double> phiS;
        std::vector<double> phitild;
        std::vector<double> phihat;
        for (k = 0; k != lT; k++)
        {
            std::default_random_engine generator;
            phiS.emplace_back(std::normal_distribution<double>{0, 1}(generator));
            phitild.emplace_back(std::normal_distribution<double>{0, 1}(generator));
            phihat.emplace_back(std::normal_distribution<double>{0, 1}(generator));
        }
        std::vector<double> S;
        std::vector<double> sigma;
        std::vector<double> r;
        std::vector<double> intr;
        S[0] = S0;
        sigma[0] = sigma0;
        r[0] = r0;
        intr[0] = r[0];
        for (k = 0; k != (lT - 1); k++)
        {
            r[k + 1] = (nu - etha * r[k]) * deltaT + sqrt(alpha * r[k]) * (rhoSr * sqrt(deltaT) * phiS[k] + sqrt(1 - rhoSr * rhoSr) * sqrt(deltaT) * phihat[k]) + r[k];
            sigma[k + 1] = kappa * (theta - sigma[k]) * deltaT + lamb * sqrt(sigma[k]) * (rhoSsigma * sqrt(deltaT) * phiS[k] + sqrt(1 - rhoSsigma * rhoSsigma) * sqrt(deltaT) * phitild[k]) + sigma[k];
            S[k + 1] = S[k] * (r[k] * deltaT + sqrt(sigma[k] * deltaT) * phiS[k]) + S[k];
            intr[k + 1] = intr[k] + r[k + 1];
        }
        paths.insert({ {"Asset",S},{"Volatility",sigma},{"Interest rate",r},{"Cum. interest rate",intr} });
    }
    return paths;
}

std::vector<double> alphaBetaTRS(int TRSpaydate, int prevTRSpaydate, typePaths paths, int nbSimul)
{
    auto alphak = 0.0;
    auto betak = 0.0;
    for (i = 0; i != nbSimul; i++)
    {
        auto res1 = (paths[i]['Asset'][TRSpaydate] - paths[i]['Asset'][prevTRSpaydate]) / (exp(paths[i]['Cum. interest rate'][TRSpaydate]));
        auto res2 = (paths[i]['Asset'][prevTRSpaydate]) / (exp(paths[i]['Cum. interest rate'][TRSpaydate]));
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

std::vector<double> TRS(std::vector<double> TRSpaymentSchedule, double sTRS, typePaths paths, int nbSimul)
    {
    auto res = 0.0;
    auto alpha = 0.0;
    auto beta = 0.0;
    for (auto s : TRSpaymentSchedule)
    {
        auto t = s--;
        alpha += alphaBetaTRS(*s, *t, paths, nbSimul)[0];
        beta += alphaBetaTRS(*s, *t, paths, nbSimul)[1] * (*s - *t);
        resTemp = alpha - sTRS * beta;
        impliedRate = alpha / beta;
    }
    std::vector<double> res;
    res.emplace_back(resTemp);
    res.emplace_back(impliedRate);
    return res;
    }

int main()
{
    auto paths = pathsGeneration(nbSimul, lT, S0, sigma0, r0, kappa, theta, lamb, nu, etha, alpha, rhoSsigma, rhoSr);
    std::cout << "TRS market value is" << TRS(TRSpaymentSchedule, sTRS, paths, nbSimul)[0] << std::endl;
    std::cout << "TRS implied rate is" << TRS(TRSpaymentSchedule, sTRS, paths, nbSimul)[1] << std::endl;
}

