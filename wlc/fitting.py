import wlc.models

import lmfit
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt

plt.style.use('ggplot')

class WormLikeChain:
    def __init__(self, model : str = "odijk") -> None:
        """Estimating the main parameters of a Worm-Like Chain Molecule from Force-Extension Measurements.
        
        Parameters
        ----------
        method : str
            Worm-like chain model for fitting. Options are: 'WLC', 'bouchiat', and 'odijk'.
        """
        # Store model functions
        if model == "WLC":
            self.func = wlc.models.WLC
        elif model == "bouchiat":
            self.func = wlc.models.bouchiat
        elif model == "odijk":
            self.func = wlc.models.odijk
        elif model == "eWLC":
            self.func = wlc.models.eWLC
        else:
            raise ValueError("Unknown fitting model. Available models are: 'WLC', 'bouchiat', and 'odijk'.")
        
        # Create fitting model
        self.fmodel = lmfit.Model(self.func)
        self.model = model

    def __repr__(self):
        """Set model representation."""
        # Recover the employed method
        return self.func.__doc__
    
    def compile(self, params : dict) -> None:
        """Compile fitting model based on parameters dictionary.
        
        Parameters
        ----------
        params : dict
            Dictionary with initial values and bounds for model parameters."""
        
        # Store general parameters
        KB = 0.013806 # Boltzmann constant in pN*nm*K-1
        if ("kBT" in params) & ("T" not in params): 
            self.fmodel.set_param_hint('kBT', value=params["kBT"])
        elif ("kBT" not in params) & ("T" in params):
            self.fmodel.set_param_hint('kBT', value=KB*(273.15 + params["T"]))
        else:
            if params["kBT"] == KB*(273.15 + params["T"]):
                self.fmodel.set_param_hint('kBT', value=params["kBT"])
            else:
                raise ValueError("Please provide either T or kBT. They are contradictory")
        
        # Store DNA parameters
        self.fmodel.set_param_hint('Lc', value=params['Lc'], min=params['Lc_lower'], max=params['Lc_upper'])
        self.fmodel.set_param_hint('Lp', value=params['Lp'], min=params['Lp_lower'], max=params['Lp_upper'])
        if (self.model == "odijk") or (self.model == "eWLC"):
            self.fmodel.set_param_hint('S', value=params['S'], min=params['S_lower'], max=params['S_upper'])

        # Create parameters
        self.fparams = self.fmodel.make_params()

        # keep kBT constant
        self.fparams['kBT'].vary = False
    
    def fit(self, data : tuple, min_delta : float, max_iters : int, filename : str = None, method : str = "leastsq" , verbose : bool = False) -> pd.DataFrame:
        """Worm-Like Chain Molecule fitting from Force-Extension Measurements.
        
        Parameters
        ----------
        data : tuple
            Observations of distance [um] and force [pN], respectively.
        min_delta : float 
            Stop fitting when increment in Lp is lower than min_delta. Units: [nm]
        max_iters : int
            Stop fitting when reach max_iters
        method : str
            Name of minimization method to use (default is "leastsq").
        filename : str
            Name of the file conatining the data. It corresponds to plot title.
        verbose : bool
            If True, a progress bar shows the progress

        Outputs
        -------
        df_file : pd.DataFrame
            Optimal results after fitting.
        """
        # Save data
        d, F = data

        # Empty dataframe to store results.
        df_full = pd.DataFrame(columns =['Lc[nm]', 'Lp[nm]','S[pN]', 'Chisqr', 'filename', 'iter', 'dLp[nm]'])

        # Stopping parameters
        dLp = np.infty
        i = 0

        # Start fitting
        with tqdm(total=max_iters, disable = not verbose) as pbar:
            while (dLp > min_delta) and (i < max_iters):
                # Fit the model
                if self.model == "odijk":
                    self.result = self.fmodel.fit(d, self.fparams, F=F, method=method)
                    Lc, Lp, S, Chisqr = self.result.params['Lc'].value, self.result.params['Lp'].value, self.result.params['S'].value, self.result.chisqr
                elif self.model == "eWLC":
                    self.result = self.fmodel.fit(F, self.fparams, fdata=data, method=method)
                    Lc, Lp, S, Chisqr = self.result.params['Lc'].value, self.result.params['Lp'].value, self.result.params['S'].value, self.result.chisqr                   
                else:
                    self.result = self.fmodel.fit(F, self.fparams, d=d, method=method)
                    Lc, Lp, S, Chisqr = self.result.params['Lc'].value, self.result.params['Lp'].value, None, self.result.chisqr

                # Check convercence
                if i > 0:
                    dLp = np.abs(Lp-df_full.loc[i-1, 'Lp[nm]'])

                df_res = pd.DataFrame({'Lc[nm]' : [Lc],
                                        'Lp[nm]' : [Lp],
                                        'S[pN]' : [S],
                                        'Chisqr' : [Chisqr],
                                        'filename' : [filename],
                                        'iter' : [i + 1],
                                        'dLp[nm]' : [dLp]})
                
                if (Lc <= self.fparams["Lc"].min) | (Lc >= self.fparams["Lc"].max):
                    print(f'Fitted Lc is out of bounds in iter {i}, data filtered out')
                elif (Lp <= self.fparams["Lp"].min) | (Lp >= self.fparams["Lp"].max):
                    print('Fitted Lp is out of bounds, data filtered out')
                else:
                    # print('Final Lc, Lp, S values within filtering bounds')
                    df_full = pd.concat((df_full, df_res), ignore_index=True)

                # Update Lp and iter
                self.fparams["Lp"].set(Lp)
                i += 1

                # Update the progress bar
                pbar.update(1)

            # Complete the progress bar
            pbar.update(max_iters-i)

            # Save fitted values per file
            df_file = pd.DataFrame([df_full[['Lc[nm]', 'Lp[nm]', 'S[pN]', 'iter', 'filename']].values[-1]], columns=['opt_Lc[nm]', 'opt_Lp[nm]', 'opt_S[pN]', 'nFittings', 'filename'])

            # Save fitting progress as attribute
            self.df_full = df_full
                
            return df_file
        
    def plot(self, data : tuple, filename : str = None):
        """Plot fitting results compare to observations.
        
        Parameters
        ----------
        data : tuple
            Observations of distance [um] and force [pN], respectively.
        """
        # Save data
        d, F = data

        # Create figure
        fig, ax = plt.subplots(figsize=(10, 7.5), dpi=300)

        if filename is not None:
            fig.suptitle('FD curve fitting (%s)' %filename)

        # Set labels
        ax.set_xlabel(r"distance [$\mu m$]")
        ax.set_ylabel(r"force [pN]")

        # Plot observations
        ax.scatter(d, F, s=2.5, c='black', label="data")

        # Plot best fit
        if self.model == "odijk":
            ax.plot(self.result.best_fit, F, c = 'red', lw = 2.5, label=self.model)
        else:
            # self.result = self.fmodel.fit(F, self.fparams, d=d)
            ax.plot(d, self.result.best_fit, c = 'red', lw = 2.5, label=self.model)
        plt.legend()
        plt.show()

    def stats(self):
        """Show fitting results."""
        print(self.result.fit_report())
    
    def plot_residuals(self, data : tuple):
        """Show residuals as a function of force or distance.
        
        Parameters
        ----------
        data : tuple
            Observations of distance [um] and force [pN], respectively.
        """
        
        # Save data
        d, F = data

        # Determine nice limits by hand
        def bins(x, y):
            binwidth = 0.25
            xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
            lim = (int(xymax/binwidth) + 1) * binwidth
            return np.arange(-lim, lim + binwidth, binwidth)

        # Create figure
        fig, axs = plt.subplots(1, 2, figsize=(7.5, 5), width_ratios=(4, 1), sharey = True, dpi=300)

        # Set labels
        if self.model == "odijk":
            axs[0].set_ylabel(r"residuals [$\mu m$]")
            axs[0].set_xlabel(r"force [pN]")
            axs[0].scatter(F, self.result.residual, s=5)
        else:
            axs[0].set_xlabel(r"distance [$\mu m$]")
            axs[0].set_ylabel(r"residuals [pN]")
            axs[0].scatter(d, self.result.residual, s=5)
        
        axs[1].hist(self.result.residual, orientation="horizontal")
        axs[1].set_xlabel("counts")
        plt.show()