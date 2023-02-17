Paul Wilmott’s program of study

I – Introduction
In his famous book http://iman.sn/bibliotek/livres/filieres/banque-finance-assurance/pdfs/paul-wilmott-on-quantitative-finance.pdf, Paul Wilmott highlighted several program of studies about numerical methods (especially finite difference and monte-carlo simulation). The purpose of this git-hub stream will be to follow this program for MC methods and show numerical results for Securities finance products, mostly for derivatives (Total return swaps)

II – Programs of study
Programs of study are following ones:

- Vanilla TRS  on a single equity
- longstaff-schwartz algorithm for cancellable TRS
- Path-dependent embedded option on a single equity : performance/financing barrier, average performance price/strike, lookback
- TRS on many stocks (Basket return swap) : Price a multi-asset contract by simulating correlated random walks
- Interest rate derivatives, spot rate model
- TRS with floating financing rate, HJM model
- TRS with floating financing rate, LMM model
- Bonus : TRS with floating financing rate, Generalized RFR MM model

III – Disclaimers
-	For each stream, a summary documentation will be attached (merely mathematical oriented)
-	Targets are Python and C++ codes, and especially numerical results
-	Sources will be documented (of course Paul Wilmott’s book but also all other sources)
-	Python and C++ codes will be assumed to be improved, but I’ll try to do my best in order to be in line with modern python and C++ concepts, without “old way” programming
