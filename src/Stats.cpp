#include "Stats.h"
#define SQ2MPI 2.506628275

double Stats::LnLogNormal(double x, double xhat, double kappa)
{
  const double lnKappa = log(kappa);
  return -log(SQ2MPI * lnKappa * x) - pow(log(x / xhat) / lnKappa, 2.0) / 2.0;
}

double Stats::LnPoisson(double n, double nexpected)
{
  return log(nexpected) * n - nexpected - lgamma(n + 1);
}
double Stats::LogLikelihood(LogLikelihoodParams *par)
{
  double logLikelihood = 0.0;
  for (size_t i = 0; i < 6; i++)
  {
    logLikelihood += LnPoisson(par->data[i], par->CRi[i] * par->transferFac);
    logLikelihood += LnLogNormal(par->CRi[i], par->CRi_expected[i], 1.0 + par->CRiErr_expected[i] / par->CRi_expected[i]);
  }
  logLikelihood += LnLogNormal(par->transferFac, par->transferFac_expected, 1.0 + par->transferFacErr_expected / par->transferFac_expected);
  return logLikelihood;
}

double Stats::minFunc(const gsl_vector *v, void *params)
{
  LogLikelihoodParams *par = (LogLikelihoodParams *)params;
  for (size_t i = 0; i < 6; i++)
  {
    par->CRi[i] = gsl_vector_get(v, i);
  }
  par->transferFac = gsl_vector_get(v, 6);
  if (par->mode == 0)
    par->signalStrength = gsl_vector_get(v, 7);

  return -Stats::LogLikelihood(par);
}

double Stats::LogLikelihoodMax(LogLikelihoodParams &par, double mu, int verbose)
{
  par.signalStrength = mu;
  par.mode = 1;
  const gsl_multimin_fminimizer_type *T =
      gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  /* Starting point */
  x = gsl_vector_alloc(7);
  for (size_t i = 0; i < 6; i++)
  {
    gsl_vector_set(x, i, par.CRi_expected[i]);
  }
  gsl_vector_set(x, 6, par.transferFac_expected);

  /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc(7);
  for (size_t i = 0; i < 6; i++)
  {
    gsl_vector_set(ss, i, par.CRiErr_expected[i]);
  }
  gsl_vector_set(ss, 6, par.transferFacErr_expected);

  /* Initialize method and iterate */
  minex_func.n = 7;
  minex_func.f = Stats::minFunc;
  minex_func.params = &par;

  s = gsl_multimin_fminimizer_alloc(T, 7);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

  do
  {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);

    if (status)
      break;

    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, 1e-6);

    if (status == GSL_SUCCESS && verbose)
    {
      printf("converged to minimum at\n");
      printf("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",
             iter,
             gsl_vector_get(s->x, 0),
             gsl_vector_get(s->x, 1),
             s->fval, size);
    }
  } while (status == GSL_CONTINUE && iter < 10000);
  if (verbose)
  {
    for (size_t i = 0; i < 6; i++)
    {
      printf("Par.Nr. %zu   %E  %E\n", i, gsl_vector_get(s->x, i), par.CRi_expected[i]);
    }
    printf("Par.Nr. %d   %E  %E\n", 6, gsl_vector_get(s->x, 6), par.transferFac_expected);
  }
  double res = s->fval;

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);

  return -res;
}

double Stats::LogLikelihoodMax(LogLikelihoodParams &par, double mu_min, double mu_max, int verbose)
{
  par.mode = 0;
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  /* Starting point */
  x = gsl_vector_alloc(8);
  for (size_t i = 0; i < 6; i++)
  {
    gsl_vector_set(x, i, par.CRi_expected[i]);
  }
  gsl_vector_set(x, 6, par.transferFac_expected);
  gsl_vector_set(x, 7, (mu_min + mu_max) / 2.0);

  /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc(8);
  for (size_t i = 0; i < 6; i++)
  {
    gsl_vector_set(ss, i, par.CRiErr_expected[i]);
  }
  gsl_vector_set(ss, 6, par.transferFacErr_expected);
  gsl_vector_set(ss, 7, (mu_min + mu_max) / 20.0);

  /* Initialize method and iterate */
  minex_func.n = 8;
  minex_func.f = Stats::minFunc;
  minex_func.params = &par;

  s = gsl_multimin_fminimizer_alloc(T, 8);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

  do
  {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    if (status)
      break;

    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, 1e-6);

    if (status == GSL_SUCCESS && verbose)
    {
      printf("converged to minimum at\n");
      printf("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",
             iter,
             gsl_vector_get(s->x, 0),
             gsl_vector_get(s->x, 1),
             s->fval, size);
    }
  } while (status == GSL_CONTINUE && iter < 10000);
  if (verbose)
  {
    for (size_t i = 0; i < 6; i++)
    {
      printf("Par.Nr. %zu   %E  %E\n", i, gsl_vector_get(s->x, i), par.CRi_expected[i]);
    }
    printf("Par.Nr. %d   %E  %E\n", 6, gsl_vector_get(s->x, 6), par.transferFac_expected);
    printf("Signalstrength  %E\n", gsl_vector_get(s->x, 7));
  }
  double res = s->fval;

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
  if (par.signalStrength < mu_min)
    return Stats::LogLikelihoodMax(par, mu_min, verbose);
  if (par.signalStrength > mu_max)
    return Stats::LogLikelihoodMax(par, mu_min, verbose);
  return -res;
}

void Stats::ReadParameters(LogLikelihoodParams &par)
{
  TFile f("FinalSamples.root", "read");
  f.Print();
  TVectorD *vSignalBins = (TVectorD *)f.Get("SigBins");
  TVectorD *vSigErr = (TVectorD *)f.Get("SigErr");
  TVectorD *vNCR_i = (TVectorD *)f.Get("NCR_i");
  TVectorD *vDataBins = (TVectorD *)f.Get("DataBins");
  TVectorD *vTransferFac = (TVectorD *)f.Get("TransferFac");
  vector<double> data, signal, CRi, CRiErr;
  for (size_t i = 0; i < 6; i++)
  {
    data.push_back((*vDataBins)[i]);
    signal.push_back((*vSignalBins)[i]);
    CRi.push_back((*vNCR_i)[i]);
    CRiErr.push_back(sqrt((*vNCR_i)[i]));
  }
  double transferFac, transferFacErr;
  transferFac = (*vTransferFac)[0];
  transferFacErr = (*vTransferFac)[1];
  par.CRi_expected = CRi;
  par.CRiErr_expected = CRiErr;
  par.transferFac_expected = transferFac;
  par.transferFacErr_expected = transferFacErr;
  par.CRi = vector<double>(6);
  par.data = data;
  par.signal = signal;
  par.signalStrength = 1;
  par.mode = 0;
}

void Stats::Test()
{
  double mu_test = 100.1;
  size_t Asimov_bg = 100;
  size_t Asimov_sig = 100;

  LogLikelihoodParams par, par_obs_0, par_obs_mu, par_pseudo;
  Stats::ReadParameters(par);
  for (size_t i = 0; i < 6; i++)
  {
    cout << par.data[i] << "  " << par.CRi_expected[i]*par.transferFac_expected<< endl;
  }
  
  par_obs_mu = par;
  par_obs_0 = par;
  double q_obs_mu = -2.0 * (Stats::LogLikelihoodMax(par_obs_mu, mu_test, 0) - Stats::LogLikelihoodMax(par, 0.0, mu_test, 0));
  Stats::LogLikelihoodMax(par_obs_0, 0.0, 0);
  std::random_device rd;
  std::mt19937 gen(rd());
  vector<poisson_distribution<int>> bgDist;
  vector<poisson_distribution<int>> sigDist;
  for (size_t i = 0; i < 6; i++)
  {
    bgDist.emplace_back(par_obs_0.CRi[i] * par_obs_0.transferFac);
    sigDist.emplace_back(par_obs_mu.CRi[i] * par_obs_mu.transferFac + mu_test * par_obs_mu.signal[i]);
  }
  int totalBG = 0;
  int totalSig = 0;
  int acceptedBG = 0;
  int acceptedSig = 0;
  double q_mu = 0.0;
  for (size_t i = 0; i < Asimov_bg; i++)
  {
    par_pseudo = par;
    for (size_t j = 0; j < 6; j++)
    {
      par_pseudo.data[j] = bgDist[j](gen);
    }
    q_mu = -2.0 * (Stats::LogLikelihoodMax(par_pseudo, mu_test, 0) - Stats::LogLikelihoodMax(par_pseudo, 0.0, mu_test, 0));
    totalBG++;
    //cout<<q_mu<<endl;
    if (q_mu > q_obs_mu)
      acceptedBG++;
  }

  for (size_t i = 0; i < Asimov_sig; i++)
  {
    par_pseudo = par;
    for (size_t j = 0; j < 6; j++)
    {
      par_pseudo.data[j] = sigDist[j](gen);
    }
    try
    {
      q_mu = -2.0 * (Stats::LogLikelihoodMax(par_pseudo, mu_test, 0) - Stats::LogLikelihoodMax(par_pseudo, 0.0, mu_test, 0));
    }
    catch (const std::exception &e)
    {
      std::cerr << e.what() << '\n';
      continue;
    }

    totalSig++;
    if (q_mu > q_obs_mu)
      acceptedSig++;
  }
  cout << q_obs_mu<< endl;
  cout << "acceptedSig " << acceptedSig << " Asimov Sig" << Asimov_sig << endl;
  cout << "acceptedBG  " << acceptedBG << " Asimov BG " << Asimov_bg << endl;
  cout << 1.0 * acceptedSig * Asimov_bg / (acceptedBG * Asimov_sig) << endl;
}