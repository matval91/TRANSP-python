    def write_NB_out(self, timeslice):
        """ write NB to ascii

        Method writing a timeslice of a simulation
        rho, power transferred to electrons and ions, the fast ion density and pressure and the current drive

        Parameters:
            None
        Attributes:
            None    
        Note:
        """
        ind = np.argmin(np.abs(self.t-timeslice))
        data = [self.rho, self.nb_FKP_vars['pe'][ind,:], self.nb_FKP_vars['pi'][ind,:],\
                self.nb_FKP_vars['n'][ind,:], self.nb_FKP_vars['press_EP'][ind,:],\
                self.nb_FKP_vars['nbcd'][ind,:]]
        header = 'rho tor\t  Pe[W/m3]\t Pi[W/m3]\t n[1/m3]\t press[Pascal]\t nbcd[A/m2]'
        out_fname='output_'+str(timeslice*100.)+'.dat'
        np.savetxt(out_fname, np.transpose(data), fmt='%.5e', delimiter='\t', header=header)

    def write_prof_in(self, timeslice):
        """ write profiles to ascii

        Method writing a timeslice of a simulation
        rho, ne, ni, Te, Ti, Zeff

        Parameters:
            None
        Attributes:
            None    
        Note:
        """
        ind = np.argmin(np.abs(self.t-timeslice))
        data = [self.rho, self.kin_vars['ne'][ind,:], self.kin_vars['ni'][ind,:], self.kin_vars['te'][ind,:], self.kin_vars['ti'][ind,:]]
        header = 'rho tor\t  ne[1/m3]\t ni[1/m3]\t Te[eV]\t\t Ti[eV]'
        out_fname='input_'+str(timeslice*100.)+'.dat'
        np.savetxt(out_fname, np.transpose(data), fmt='%.5e', delimiter='\t', header=header)
