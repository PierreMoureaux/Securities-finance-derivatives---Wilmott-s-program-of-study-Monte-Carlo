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
auto r = 0.05;

//TRS datas
auto T = 1.0;
std::vector<double> TRSpaymentSchedule = { 0,25, 50, 75, 99 };
auto sTRS = 0.01;

//Model datas
auto kappa = 0.1;
auto theta = 0.1;
auto gamma = 0.2;
auto rhoSsigma = 0.1;
auto lamb = 20;
auto muJ = 0.1;
auto sigmaJ = 0.1;
auto deltaT = 0.01;
auto lT = static_cast<int>(T / deltaT);
auto nbSimul = 1000;

typePaths pathsGeneration(int nbSimul, int lT, double S0, double sigma0, double r, double kappa, double theta, double gamma, 
    double lamb, double muJ, double sigmaJ, double rhoSsigma)
{
    typePaths paths;
    for (auto i = 0; i != nbSimul; i++)
    {
        std::vector<double> phiS;
        std::vector<double> phitild;
        std::vector<double> phiJ;
        std::vector<double> P;
        std::default_random_engine generator;
        for (auto k = 0; k != lT; k++)
        {
            phiS.emplace_back(std::normal_distribution<double>{0, 1}(generator));
            phitild.emplace_back(std::normal_distribution<double>{0, 1}(generator));
            phiJ.emplace_back(std::normal_distribution<double>{muJ, sigmaJ}(generator));
            P.emplace_back(std::poisson_distribution<int>{deltaT*lamb}(generator));
        }
        std::vector<double> S(lT);
        std::vector<double> sigma(lT);
        std::vector<double> intr(lT);
        S[0] = S0;
        sigma[0] = sigma0;
        intr[0] = r;
        auto meanJRN = exp(muJ + 0.5 * (sigmaJ*sigmaJ)) - 1;
        if (2 * kappa * theta <= (gamma*gamma))
        {
            kappa = (gamma*gamma) / (2 * theta) + 5;
        }
        for (auto k = 0; k != (lT - 1); k++)
        {
            sigma[k + 1] = kappa * (theta - sigma[k]) * deltaT + gamma * sqrt(sigma[k]) * (rhoSsigma * sqrt(deltaT) * phiS[k] + sqrt(1 - rhoSsigma * rhoSsigma) * sqrt(deltaT) * phitild[k]) + sigma[k];
            S[k + 1] = S[k] * ((r - lamb * meanJRN)* deltaT + sqrt(sigma[k] * deltaT) * phiS[k] + (exp(phiJ[k]) - 1) * P[k]) + S[k];
            intr[k + 1] = intr[k] + r;
        }
        std::unordered_map<std::string, std::vector<double>> interRes = { {"Asset",S},{"Volatility",sigma},{"Cum. interest rate",intr} };
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
    auto paths = pathsGeneration(nbSimul, lT, S0, sigma0, r, kappa, theta, gamma,lamb, muJ, sigmaJ, rhoSsigma);
    auto TRSres = TRS(TRSpaymentSchedule, sTRS, paths, nbSimul);
    std::cout << "TRS market value is" << TRSres[0] << std::endl;
    std::cout << "TRS implied rate is" << TRSres[1] << std::endl;
}
