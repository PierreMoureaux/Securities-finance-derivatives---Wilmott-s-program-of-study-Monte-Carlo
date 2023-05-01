#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <math.h>
#include <any>
#include <functional>
#include <numeric>
#include <algorithm>

typedef std::vector<std::unordered_map<std::string, std::vector<double>>> typePaths;

//Market datas
auto r0 = 0.05;
auto lambda0 = 0.1;
auto gamma0 = 0.05;

//TRS datas
auto T = 1;
std::vector<double> TRSpaymentSchedule = { 0,25, 50, 75, 99 };
auto sTRS = 0.02;

//Bond datas
auto TP = 2;
std::vector<double> couponSchedule = { 13, 33, 63 };
auto c = 0.03;
auto R = 0.4;

//Model datas
auto nu = 0.1;
auto etha = 0.2;
auto vola = 0.1;
auto alpha = 0.1;
auto beta = 0.2;
auto sigma = 0.1;
auto omega = 0.1;
auto deltaT = 0.01;
auto lT = static_cast<int>(T / deltaT);
auto lTP = static_cast<int>(TP / deltaT);
auto nbSimul = 100;
auto nbSimulP = 100;

double BondPrice(double nbSimulP, int tk, const std::vector<double>& couponSchedule, double c, int lTP, double r0, double lambda0,
    double gamma0, double nu, double etha, double vola, double alpha, double beta, double sigma, double omega, double R)
{
    auto Ptk = 0.0;
    std::vector<double> nestedPath(nbSimulP);
    for (auto i = 0; i != nbSimulP; i++)
    {
        std::vector<double> phir;
        std::vector<double> phil;
        std::vector<double> phig;
        std::default_random_engine generator;
        for (auto k = 0; k != lTP; k++)
        {
            phir.emplace_back(std::normal_distribution<double>{0, 1}(generator));
            phil.emplace_back(std::normal_distribution<double>{0, 1}(generator));
            phig.emplace_back(std::normal_distribution<double>{0, 1}(generator));
        }
        std::vector<double> r(lTP);
        std::vector<double> l(lTP);
        std::vector<double> g(lTP);
        auto c1 = 0.0;
        auto c2 = 0.0;
        auto c3 = 0.0;
        r[0] = r0;
        l[0] = lambda0;
        g[0] = gamma0;
        for (auto k = 0; k != (lTP - 1); k++)
        {
            r[k + 1] = etha * (nu - r[k]) * deltaT + vola * sqrt(r[k] * deltaT) * phir[k] + r[k];
            l[k + 1] = beta * (alpha - l[k]) * deltaT + sigma * sqrt(l[k] * deltaT) * phil[k] + l[k];
            g[k + 1] = omega * sqrt(deltaT) * phig[k] + g[k];
        }
        c2 = exp(-(std::accumulate(r.begin() + tk, r.end(), 0) + std::accumulate(l.begin() + tk, l.end(), 0) + std::accumulate(g.begin() + tk, g.end(), 0)));
        for (auto t = (tk + 1); t != (lTP - 1); t++)
        {
            auto lastDate = (lTP - 1 - t);
            if (std::find(couponSchedule.begin(), couponSchedule.end(), t) != couponSchedule.end())
            {
                c1 += exp(-(std::accumulate(r.begin() + tk, r.end() - lastDate, 0) + std::accumulate(l.begin() + tk, l.end() - lastDate, 0) + std::accumulate(g.begin() + tk, g.end() - lastDate, 0)));
            }
            c3 += l[t]*exp(-(std::accumulate(r.begin() + tk, r.end() - lastDate, 0) + std::accumulate(l.begin() + tk, l.end() - lastDate, 0) + std::accumulate(g.begin() + tk, g.end() - lastDate, 0)));
        }
        nestedPath[i] = c * c1 + c2 + (1 - R) * c3;
        Ptk += nestedPath[i];
    }
    return Ptk / nbSimulP;
}

typePaths pathsGeneration(int nbSimul, int nbSimulP, int lT, const std::vector<double>& couponSchedule, double c, int lTP, double r0, double lambda0,
    double gamma0, double nu, double etha, double vola, double alpha, double beta, double sigma, double omega, double R)
{
    typePaths paths;
    for (auto i = 0; i != nbSimul; i++)
    {
        std::default_random_engine generator;
        std::vector<double> phir2;
        for (auto k = 0; k != lT; k++)
        {
            std::default_random_engine generator;
            phir2.emplace_back(std::normal_distribution<double>{0, 1}(generator));
        }
        std::vector<double> P(lT);
        std::vector<double> r(lT);
        std::vector<double> intr(lT);
        r[0] = r0;
        intr[0] = r[0];
        for (auto k = 0; k != (lT - 1); k++)
        {
            r[k + 1] = etha * (nu - r[k]) * deltaT + vola * sqrt(r[k] * deltaT) * phir2[k] + r[k];
            intr[k + 1] = intr[k] + r[k + 1];
            P[k] = BondPrice(nbSimulP, k, couponSchedule, c, lTP, r0, lambda0, gamma0, nu, etha, vola, alpha, beta, sigma, omega, R);
        }
        std::unordered_map<std::string, std::vector<double>> interRes = {{"Asset",P}, {"Interest rate",r}, {"Cum. interest rate",intr}};
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

double cTRS(int TRScoupondate, const typePaths& paths, int nbSimul)
{
    auto ck = 0.0;
    for (auto i = 0; i != nbSimul; i++)
    {
        auto res3 = 1.0 / (exp(paths[i].at("Cum. interest rate")[TRScoupondate]));
        ck += res3;
    }
    ck /= nbSimul;
    return ck;
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

std::vector<double> TRS(const std::vector<double>& TRSpaymentSchedule, 
    const std::vector<double>& couponSchedule, double sTRS, const typePaths& paths, int nbSimul)
{
    auto alpha = 0.0;
    auto beta = 0.0;
    auto ceta = 0.0;
    for (auto s = TRSpaymentSchedule.begin(); s != TRSpaymentSchedule.end() - 1; s++)
    {
        auto payDate = getIndex(TRSpaymentSchedule, *s) + 1;
        auto prevPayDate = getIndex(TRSpaymentSchedule, *s);
        auto abRes = alphaBetaTRS(TRSpaymentSchedule[payDate], TRSpaymentSchedule[prevPayDate], paths, nbSimul);
        alpha += abRes[0];
        beta += abRes[1] * (TRSpaymentSchedule[payDate] - TRSpaymentSchedule[prevPayDate]);
    }
    for (auto s = couponSchedule.begin(); s != couponSchedule.end() - 1; s++)
    {
        auto couponDate = getIndex(couponSchedule, *s);
        auto cRes = cTRS(couponDate, paths, nbSimul);
        ceta += cRes;
    }
    auto resTemp = alpha+c*ceta - sTRS * beta;
    auto impliedRate = (alpha+c*ceta) / beta;
    std::vector<double> res;
    res.emplace_back(resTemp);
    res.emplace_back(impliedRate);
    return res;
}

int main()
{
    auto paths = pathsGeneration(nbSimul, nbSimulP, lT, couponSchedule, c, lTP, r0, lambda0, gamma0, nu, etha, vola, alpha, beta, sigma, omega, R);
    auto TRSres = TRS(TRSpaymentSchedule, couponSchedule, sTRS, paths, nbSimul);
    std::cout << "TRS market value is" << TRSres[0] << std::endl;
    std::cout << "TRS implied rate is" << TRSres[1] << std::endl;
}


