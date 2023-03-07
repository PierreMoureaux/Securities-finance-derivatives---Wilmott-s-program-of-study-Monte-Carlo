#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <math.h>
#include <any>
#include <iostream>

typedef std::vector<std::unordered_map<std::string, std::vector<double>>> typePaths;

//Market datas
auto S0 = 100.0;
auto sigma0 = 0.1;
auto r0 = 0.05;

//TRS datas
auto T = 1.0;
std::vector<double> TRSpaymentSchedule = { 0,25, 50, 75, 99 };
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
auto nbSimul = 100;

typePaths pathsGeneration(int nbSimul, int lT, double S0, double sigma0, double r0, double kappa, double theta, double lamb, double nu, double etha, double alpha, double rhoSsigma, double rhoSr)
{
    typePaths paths;
    for (auto i = 0; i != nbSimul; i++)
    {
        std::vector<double> phiS;
        std::vector<double> phitild;
        std::vector<double> phihat;
        for (auto k = 0; k != lT; k++)
        {
            std::default_random_engine generator;
            phiS.emplace_back(std::normal_distribution<double>{0, 1}(generator));
            phitild.emplace_back(std::normal_distribution<double>{0, 1}(generator));
            phihat.emplace_back(std::normal_distribution<double>{0, 1}(generator));
        }
        std::vector<double> S(lT);
        std::vector<double> sigma(lT);
        std::vector<double> r(lT);
        std::vector<double> intr(lT);
        S[0] = S0;
        sigma[0] = sigma0;
        r[0] = r0;
        intr[0] = r[0];
        for (auto k = 0; k != (lT - 1); k++)
        {
            r[k + 1] = (nu - etha * r[k]) * deltaT + sqrt(alpha * r[k]) * (rhoSr * sqrt(deltaT) * phiS[k] + sqrt(1 - rhoSr * rhoSr) * sqrt(deltaT) * phihat[k]) + r[k];
            sigma[k + 1] = kappa * (theta - sigma[k]) * deltaT + lamb * sqrt(sigma[k]) * (rhoSsigma * sqrt(deltaT) * phiS[k] + sqrt(1 - rhoSsigma * rhoSsigma) * sqrt(deltaT) * phitild[k]) + sigma[k];
            S[k + 1] = S[k] * (r[k] * deltaT + sqrt(sigma[k] * deltaT) * phiS[k]) + S[k];
            intr[k + 1] = intr[k] + r[k + 1];
        }
        typePaths interRes = { {"Asset",S},{"Volatility",sigma},{"Interest rate",r},{"Cum. interest rate",intr} };
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
        auto res1 = (paths[i].at("Asset")[TRSpaydate] - paths[i].at("Asset")[prevTRSpaydate]) / (exp(paths[i].at("Cum. interest rate")[TRSpaydate]));
        auto res2 = (paths[i].at("Asset")[prevTRSpaydate]) / (exp(paths[i].at("Cum. interest rate")[TRSpaydate]));
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
    for (auto s = TRSpaymentSchedule.begin();s!= TRSpaymentSchedule.end()-1;s++)
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
    auto paths = pathsGeneration(nbSimul, lT, S0, sigma0, r0, kappa, theta, lamb, nu, etha, alpha, rhoSsigma, rhoSr);
    auto TRSres = TRS(TRSpaymentSchedule, sTRS, paths, nbSimul);
    std::cout << "TRS market value is" << TRSres[0] << std::endl;
    std::cout << "TRS implied rate is" << TRSres[1] << std::endl;
}