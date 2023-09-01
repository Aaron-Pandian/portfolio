# Halcyon Stage Prop Calculator

import csv
import json
import numpy as np
import CoolProp.CoolProp as CP
from scipy.optimize import newton
import scipy.signal
import sys
import yaml
import time

import pandas as pd
import os
import nptdms
import matplotlib.pyplot as pl

mdot_LOX_RP1_pers = []


class Halcyon:
    def __init__(self):
        # Constant Inputs ----------------------------------------------------------
        self.R = 8.314441
        self.pi = 3.1415
        self.g0 = 9.81

        # Press Fluid
        self.fluidPress = 'Helium'
        self.maxImpulse = 200000  # lbf-s

        # Combustion Efficiency Inputs
        self.Cstar_eff = 0.95

        # Time Loop Variables
        self.dT = 0.1  # seconds per index
        self.Pconfig = " "

        # Unit Conversions
        self.L2m3 = 0.001
        self.m32L = 1 / 0.001
        self.pa2psi = 1 / 6895
        self.psi2pa = 6895
        self.in2m = 0.0254
        self.m2in = 39.3701
        self.N2lbf = 0.224809
        self.lbf2N = 1 / 0.224809

        # Rocket heights from engine cap plane in m measured to bottom of tanks, from P&ID
        self.x0_RP1_tank = 0.76 - 0.1515  # subtracting dome height from EB skirt
        self.dx_RP1_tank = 1.455  # total potential fluid height in RP tank
        self.dx_OF = 0.61 - 0.1515*2  # subtracting fuel dome and lox dome from skirt height
        self.dx_LOX_tank = 2.435  # total potential fluid height in LOX tank
        self.dx_PO = 0.61 - 0.1515*2  # subtracting LOX dome and press dome from skirt height
        self.x0_LOX_tank = self.x0_RP1_tank + self.dx_RP1_tank + self.dx_OF
        self.x0_GHE_tank = self.x0_LOX_tank + self.dx_LOX_tank + self.dx_PO

        # Tank Volumes (assumes cylindrical tank with hemisphere domes), locked in
        self.V_LOX_Tank = 249 * self.L2m3
        self.r_LOX_Tank = 0.18266  # m
        self.V_RP1_Tank = 141.2 * self.L2m3
        self.r_RP1_Tank = 0.17971  # m
        self.V_GHE_Tank = 167 * self.L2m3
        self.r_GHE_Tank = 0.18101  # m
        self.V_RP1_Runline = ((1 - 2 * 0.065) * self.in2m)**2 / \
            4 * self.x0_RP1_tank  # aluminum 1" tubing
        self.V_LOX_Runline = (
            (1 - 2 * 0.049) * self.in2m)**2 / 4 * self.x0_LOX_tank

        # Press / RCS CdA Inputs, Sizing of press and RCS orifices, adjust to maintian proper press authority
        self.Cd_press_LOX = 0.61
        self.A_press_LOX = (0.11 * self.in2m)**2 / 4 * self.pi

        self.Cd_press_RP1 = 0.61
        self.A_press_RP1 = (0.065 * self.in2m)**2 / 4 * self.pi

        self.Cd_press_RCS = 0.61
        self.A_press_RCS = (0.125 * self.in2m)**2 / 4 * self.pi

        # Main Propline CdA Inputs, based on engine constants, should not be varied other than movmean to damped instabilities
        self.mdot_RP1_const = 1.8
        self.mdot_LOX_const = 1.8 * 2.26
        self.MEOP_Pc = 500
        self.A_thrt = (1.23 * self.in2m)**2 * self.pi
        self.A_exit = (4.265 * self.in2m)**2 * self.pi
        self.expan_ratio = self.A_exit / self.A_thrt
        self.n_movmean = 2

        # Rocket Drymass
        # from mass budget https://utexas.sharepoint.com/:x:/s/TREL/EYfVW7DDqJ9PtFU2T6Khx-gBwBiNDR0UVot4f7J14v--fQ?e=OtWTz2
        self.kg_drymass = 286.610

        # CEA Lookup Tables for Kerolox Combustion
        self.klookup = np.loadtxt(open("CEA/klookup.csv", "rb"), delimiter=",")
        self.Rlookup = np.loadtxt(open("CEA/Rlookup.csv", "rb"), delimiter=",")
        self.Tlookup = np.loadtxt(open("CEA/Tlookup.csv", "rb"), delimiter=",")

        # Atmospheric Properties Lookup CSV http://www.braeunig.us/space/atmos.htm
        self.atm_z_lookup = np.loadtxt(open("ATM/z.csv", "rb"), delimiter=",")
        self.atm_rho_lookup = np.loadtxt(
            open("ATM/rho.csv", "rb"), delimiter=",")
        self.atm_P_lookup = np.loadtxt(open("ATM/P.csv", "rb"), delimiter=",")

        # Custom Defined Fuel Properties
        self.fluid_data_RP1 = {
            'CAS': '8008-20-6',  # CAS Number from random kerosene page on google
            'Tc': 673.15,  # Critical Temperature from Jet A https://apps.dtic.mil/sti/pdfs/AD1093317.pdf
            'Tc_units': 'K',
            'acentric': 0.514,  # Acentric Factor from JP-5 https://apps.dtic.mil/sti/pdfs/AD1093317.pdf
            'aliases': ['KER'],
            'molemass': 0.170,  # Molar Mass
            'molemass_units': 'kg/mol',
            'name': 'kerosene',
            'pc': 2381138.0,  # Critical Pressure from Jet A https://apps.dtic.mil/sti/pdfs/AD1093317.pdf
            'pc_units': 'Pa'
        }

        # Prop Load Inits
        self.rho0_LOX = 0
        self.rho0_RP1 = 0
        self.V0_LOX = 0
        self.V0_LOX_ull = 0
        self.V0_RP1 = 0
        self.V0_RP1_ull = 0
        self.V_LOX_dome = 0
        self.V_LOX_cyl = 0
        self.V_RP1_dome = 0
        self.V_RP1_cyl = 0
        self.x0_LOX = 0
        self.x0_RP1 = 0
        self.kg0_rkt = 0

        # Varying tank property inits
        self.P_GHE = [0] * 1
        self.T_GHE = [0] * 1
        self.P_RP1 = [0] * 1
        self.P_LOX = [0] * 1
        self.kg_GHE = [0] * 1
        self.time = [0] * 1

        # Variable Inputs ----------------------------------------------------------

        # Prop Load Config (temps in kelvin, pressure in psi, mass in kg)
        self.T0_LOX = 100
        self.kg0_LOX = 221
        self.MEOP_LOX = 740

        self.T0_RP1 = 300
        self.kg0_RP1 = 96
        self.MEOP_RP1 = 855

        # Press fill configs
        self.T0_GHE = 300
        self.MEOP_GHE = 6000

        # Control System Values ----------------------------------------------------------

        # System MAWPs in psi
        self.MAWP_LOX = 740
        self.MAWP_RP1 = 1313
        self.MAWP_ENG = 550

        # Bang Bang Targets
        self.Pmin_LOX = self.MEOP_LOX - 25
        self.Pmax_LOX = self.MEOP_LOX + 25
        self.Pflag_LOX = 0  # 1 = pressing in window, 0 = no press

        self.Pmin_RP1 = self.MEOP_RP1 - 25
        self.Pmax_RP1 = self.MEOP_RP1 + 25
        self.Pflag_RP1 = 0  # 1 = pressing in window, 0 = no press

        self.Rollmax = 0
        self.Rollmin = 0
        self.Rollflag = 0  # 1 = firing RCS, 0 = no RCS firing

        # Hold Duration Post Prepress
        self.prepress_hold = 15  # in seconds

        # Tank Prop margin
        self.n_shutoff = 0.02  # MECO when any prop tank reaches <2% total volume

    # Helper Functions  ----------------------------------------------------------

    # Get density in kg/m^3 for pressure in psi, Temp in kelvin
    def get_rho(self, T_k, P_psi, Fluid):
        P_pa = P_psi * self.psi2pa
        rho = CP.PropsSI('D', 'T', T_k, 'P', P_pa, Fluid)
        return rho

    # kg cost per change in pressure for press (cf is collapse factor) - Redlich-Kwong
    def get_kg4dP(self, vol, P, T_k, cf, gas):

        Tc = CP.PropsSI("Tcrit", gas)  # critical pressure in K
        pc = CP.PropsSI("pcrit", gas) * 1000  # Critical Pressure in Pa
        P_pa = P * self.psi2pa

        a = 0.42748 * self.R ** 2 * Tc ** (5 / 2) / pc
        b = 0.08664 * self.R * Tc / pc

        a1 = P_pa
        a2 = -(self.R * T_k)
        a3 = P_pa*b**2 - self.R * T_k * b + a / np.sqrt(T_k)
        a4 = -a * b / (np.sqrt(T_k))

        coeff = [a1, a2, a3, a4]
        Vm = np.roots(coeff)
        Vm = Vm[0]
        n = vol/Vm
        mm = CP.PropsSI("molar_mass", gas)

        kg = n * mm * cf
        return kg

    # dP change for increase in ullage volume Redlich-Kwong
    def get_P4V(self, V, kg, T, cf, gas):

        Tc = CP.PropsSI("Tcrit", gas)  # critical pressure in K
        pc = CP.PropsSI("pcrit", gas) * 1000  # Critical Pressure in Pa

        a = 0.42748 * self.R ** 2 * Tc ** (5 / 2) / pc
        b = 0.08664 * self.R * Tc / pc

        mm = CP.PropsSI("molar_mass", gas)
        n = kg / (mm * cf)
        Vm = V / n

        p = self.R * T / (Vm - b) - a / (np.sqrt(T) * Vm * (Vm + b))
        p_psi = p * self.pa2psi
        return p_psi

    def get_P4kg(self, V, kg, T, cf, gas):
        # https://en.wikipedia.org/wiki/Real_gas  Redlich-Kwong

        Tc = CP.PropsSI("Tcrit", gas)  # critical pressure in K
        pc = CP.PropsSI("pcrit", gas) * 1000  # Critical Pressure in Pa

        a = 0.42748 * self.R**2 * Tc ** (5/2) / pc
        b = 0.08664 * self.R * Tc / pc

        mm = CP.PropsSI("molar_mass", gas)
        n = kg / (mm * cf)
        Vm = V / n

        p = self.R * T / (Vm - b) - a / (np.sqrt(T) * Vm * (Vm + b))
        p_psi = p * self.pa2psi
        return p_psi

    # start with lox and RP press
    def set_pressConfig(self, P_LOX, P_RP1, Roll_Hz):
        # First check flags to see if we are pressing or not from last iteration

        # Actively pressing LOX in last iteration
        if self.Pflag_LOX == 1:

            # Check to see if pressure has exceeded max bang bang, if so, kill flag
            if P_LOX > self.Pmax_LOX:
                self.Pflag_LOX = 0
            else:
                self.Pflag_LOX = self.Pflag_LOX

        # No LOX actively pressing in last iteration, look to see if any pressures fall below min bang targets
        if self.Pflag_LOX == 0:
            # Check to see if pressure has fallen below min bang bang, if so, open flag
            if P_LOX < self.Pmin_LOX:
                self.Pflag_LOX = 1
            else:
                self.Pflag_LOX = self.Pflag_LOX

        # Actively pressing RP1 in last iteration
        if self.Pflag_RP1 == 1:

            # Check to see if pressure has exceeded max bang bang, if so, kill flag
            if P_RP1 > self.Pmax_RP1:
                self.Pflag_RP1 = 0
            else:
                self.Pflag_RP1 = self.Pflag_RP1

        # No RP1 actively pressing in last iteration, look to see if any pressures fall below min bang targets
        if self.Pflag_RP1 == 0:
            # Check to see if pressure has fallen below min bang bang, if so, open flag
            if P_RP1 < self.Pmin_RP1:
                self.Pflag_RP1 = 1
            else:
                self.Pflag_RP1 = self.Pflag_RP1

        # Actively firing RCS in last iteration
        if self.Rollflag == 1:

            # Check to see if roll has dipped below min bang bang, if so, kill flag
            if Roll_Hz < self.Rollmin:
                self.Rollflag = 0
            else:
                self.Rollflag = self.Rollflag

        # No RCS actively firing in last iteration, look to see if any roll above max bang targets
        if self.Pflag_RP1 == 0:
            # Check to see if pressure has fallen below min bang bang, if so, open flag
            if Roll_Hz > self.Rollmax:
                self.Rollflag = 1
            else:
                self.Rollflag = self.Rollflag

    # Finds compressible mass flow for given fluid properties
    def get_comprMdot(self, Phigh_psi, Plow_psi, Cd, A_m2, T_k, Fluid):

        # Find Gamma (ratio of specific heats) at given press temp and pressure
        P_pa = Phigh_psi * self.psi2pa
        Cp = CP.PropsSI('C', 'T', T_k, 'P', P_pa, Fluid)
        Cv = CP.PropsSI('O', 'T', T_k, 'P', P_pa, Fluid)
        gamma = Cp / Cv

        rho = self.get_rho(T_k, Phigh_psi, Fluid)

        # Find Critical Pressure from gamma
        crit_P = Phigh_psi * ((2 / (gamma + 1)) ** (gamma / (gamma - 1)))

        # Unchoked Flow condition
        if Plow_psi > crit_P:
            mdot = Cd * A_m2 * np.sqrt(2 * Phigh_psi * self.psi2pa * rho * (gamma / (gamma - 1)) * (
                ((Plow_psi / Phigh_psi) ** (2 / gamma)) - ((Plow_psi / Phigh_psi) ** ((gamma + 1) / gamma))))

        # Choked FLow Condition
        elif Plow_psi < crit_P:
            mdot = Cd * A_m2 * np.sqrt(gamma * Phigh_psi * self.psi2pa * rho * ((2 / (gamma + 1)) **
                                                                                ((gamma + 1) / (gamma - 1))))
        return mdot

    # Kilogram shuffling of helium based on bang bang flags
    def get_pressdKg(self, P_LOX, P_RP1, P_ATM, P_GHE, T_GHE):

        # Change in kg of helium from press tank to LOX tank
        if self.Pflag_LOX == 1:
            mdot_LOX = self.get_comprMdot(
                P_GHE, P_LOX, self.Cd_press_LOX, self.A_press_LOX, T_GHE, self.fluidPress)
            dkg_LOX = mdot_LOX * self.dT
        else:
            dkg_LOX = 0

        # Change in kg of helium from press tank to RP1 tank
        if self.Pflag_RP1 == 1:
            mdot_RP1 = self.get_comprMdot(
                P_GHE, P_RP1, self.Cd_press_RP1, self.A_press_RP1, T_GHE, self.fluidPress)
            dkg_RP1 = mdot_RP1 * self.dT
        else:
            dkg_RP1 = 0

        # Change in kg of helium from press tank out rcs thrusters
        if self.Rollflag == 1:
            mdot_RCS = self.get_comprMdot(
                P_GHE, P_ATM, self.Cd_press_RCS, self.A_press_RCS, T_GHE, self.fluidPress)
            dkg_RCS = mdot_RCS * self.dT
        else:
            dkg_RCS = 0

        dkg_GHE = -1*(dkg_LOX + dkg_RP1 + dkg_RCS)

        return dkg_LOX, dkg_RP1, dkg_RCS, dkg_GHE

    # Equalize Tank Pressures from an input kg transfer

    def get_pressEq(self, V_LOX_ull, kg_LOX, V_RP1_ull, kg_RP1, T_GHE, kg_GHE):
        # Change in lox ullage pressure due to added kg
        P_LOX_next = self.get_P4kg(
            V_LOX_ull, kg_LOX, T_GHE, 2, self.fluidPress)  # Collapse factor of 2

        # New RP1 ullage pressure due to added kg
        P_RP1_next = self.get_P4kg(
            V_RP1_ull, kg_RP1, T_GHE, 1, self.fluidPress)

        # New helium pressure due to lost kgs
        P_GHE_next = self.get_P4kg(
            self.V_GHE_Tank, kg_GHE, T_GHE, 1, self.fluidPress)

        # Find New Helium Density from vol and mass
        rho_GHE = kg_GHE / self.V_GHE_Tank

        # Find new ullage temp at the lowered density and pressure
        T_GHE_next = CP.PropsSI('T', 'D', rho_GHE, 'P',
                                P_GHE_next * self.psi2pa, self.fluidPress)

        # WIP constant
        c = 50
        c2 = self.dT / c
        dT = c2 * (T_GHE_next - T_GHE)
        T_GHE_next = T_GHE + dT

        return P_LOX_next, P_RP1_next, P_GHE_next, T_GHE_next

    # Find change in ullage volume from input prop mass flow temp and pressure
    def get_dVullage(self, mdot_LOX, P_LOX, T_LOX, V_LOX_ullage, mdot_RP1, P_RP1, T_RP1, V_RP1_ullage):
        dkg_LOX = mdot_LOX * self.dT
        rho_LOX = self.get_rho(T_LOX, P_LOX, "Oxygen")
        dV_LOX = dkg_LOX / rho_LOX
        V_LOX_ullage_next = V_LOX_ullage + dV_LOX
        # caps ullage volume at tank volume
        if V_LOX_ullage_next > self.V_LOX_Tank:
            V_LOX_ullage_next = self.V_LOX_Tank
        else:
            V_LOX_ullage_next = V_LOX_ullage_next

        dkg_RP1 = mdot_RP1 * self.dT
        rho_RP1 = self.get_rho(T_RP1, P_RP1, "PR::kerosene")
        dV_RP1 = dkg_RP1 / rho_RP1
        V_RP1_ullage_next = V_RP1_ullage + dV_RP1
        # caps ullage volume at tank volume
        if V_RP1_ullage_next > self.V_RP1_Tank:
            V_RP1_ullage_next = self.V_RP1_Tank
        else:
            V_RP1_ullage_next = V_RP1_ullage_next

        return V_LOX_ullage_next, V_RP1_ullage_next

    # Finds prop mdots from Pc, tank pressure, and CdAs
    def get_propMdot(self, Pc, Ptank, Fluid):

        if Fluid == "Oxygen":
            CdA = self.CdA_LOX
            T_k = self.T0_LOX
            P_psi = (Pc + Ptank)/2
            rho = self.get_rho(T_k, P_psi, Fluid)

            # incompressible CdA equation using entire runline CdA
            mdot = CdA * np.sqrt(2*rho * (Ptank - Pc) *
                                 self.psi2pa) * (1 / self.Cstar_eff)

        elif Fluid == "PR::kerosene":
            CdA = self.CdA_RP1
            T_k = self.T0_RP1
            P_psi = (Pc + Ptank) / 2
            rho = self.get_rho(T_k, P_psi, Fluid)

            # incompressible CdA equation using entire runline CdA
            mdot = CdA * np.sqrt(2 * rho * (Ptank - Pc) *
                                 self.psi2pa) * (1 / self.Cstar_eff)
        else:
            print('ERROR: Invalid Fluid in get_propMdot!!: ', Fluid)
            mdot = 0
        return mdot

    # Gets gamma, Rgas constant, and Tflame value from OtF and Pc, models combustion from prop inputs / last state info
    def get_CEA(self, OtF, Pc):
        OtF_list = self.klookup[0, :]
        Pc_list = self.klookup[:, 0]

        jmaxP = len(Pc_list) - 1
        jmaxO = len(OtF_list) - 1

        if Pc < Pc_list[jmaxP]:
            j = 0
            while Pc > Pc_list[j]:
                self.Pc_i = j + 1
                j += 1
        else:
            self.Pc_i = 0

        if OtF < OtF_list[jmaxO]:
            j = 0
            while OtF > OtF_list[j]:
                self.OtF_i = j + 1
                j += 1
        else:
            self.OtF_i = 0

        if (self.Pc_i * self.OtF_i) == 0:
            self.Pc_i = 0
            self.OtF_i = 0

        gamma = float(self.klookup[self.OtF_i, self.Pc_i])
        gasConst = self.Rlookup[self.OtF_i, self.Pc_i]
        T_flame = self.Tlookup[self.OtF_i, self.Pc_i]

        return gamma, gasConst, T_flame

    # Finds next iteration chamber pressure from CEA inputs
    def get_Pc4CEA(self, gamma, gasConst, T_flame, mdot_tot):
        Pc = 1 / (((self.A_thrt / T_flame ** 0.5) * (gamma / gasConst) ** 0.5 *
                   (1 + (gamma - 1) / 2) ** (-1 * ((gamma + 1) / (2 * (gamma - 1))))) / mdot_tot)
        # Adds the combustion efficiency to drop chamber pressure
        Pc_psi = Pc * self.pa2psi
        return Pc_psi

    # finds thrust from CEA -> exit Mach -> exit Temp -> exit Pressure -> exit Velocity -> Thrust
    # https://www.grc.nasa.gov/www/k-12/VirtualAero/BottleRocket/airplane/isentrop.html
    def get_thrust4CEA(self, gamma, gasConst, T_flame, P_amb, Pc, mdot_tot):
        # Find exit mach, code stolen from Serenity and jerry rigged to take single value inputs
        gamma = [gamma] * 1
        gammaMe = np.array(gamma)
        gammaMe[gammaMe == 0] = 'nan'

        x = np.empty(len(gammaMe))
        x.fill(3)

        def M_exit_Calc(M_exit, gamma): return (1 / M_exit) * np.power(
            (1 + 0.5 * (gamma - 1) * M_exit ** 2) / (1 + 0.5 * (gamma - 1)),
            (gamma + 1) / (2 * (gamma - 1))) - self.expan_ratio
        M_exit = newton(M_exit_Calc, x, args=(gammaMe,), maxiter=100)
        M_exit[np.isnan(M_exit)] = 0
        M_exit = M_exit[0]
        gamma = gamma[0]

        # Find T_exit from mach, gamma, flame temp
        T_exit = T_flame * (1 + (gamma - 1) / 2 * M_exit ** 2) ** -1

        # Finds exit pressure in Pa
        P_exit_Pa = Pc * self.psi2pa * \
            (1 + (gamma - 1) / 2 * M_exit ** 2) ** -(gamma / (gamma - 1))
        P_exit = P_exit_Pa * self.pa2psi

        # Finds Exit velocity in m/s
        V_exit = M_exit * (gamma * gasConst * T_exit) ** 0.5

        thrust_N = mdot_tot * V_exit + \
            (P_exit - P_amb) * self.psi2pa * self.A_exit
        thrust = thrust_N * self.N2lbf

        Isp = thrust_N / (self.g0 * mdot_tot)

        return thrust, Isp

    # finds a moving mean of the data list inputted, n is the window of indexes it averages
    def get_movmean(self, list2, n):
        listz = [i for i in list2 if i != 0]
        if len(listz) > n:
            a = list2
            ret = np.cumsum(a, dtype=float)
            ret[n:] = ret[n:] - ret[:-n]
            list2 = ret[n - 1:] / n

            size = len(list2) + n - 1
            list2.resize(size, )

            for i in range(len(list2)-(n-1), len(list2)):
                list2[i] = list2[len(list2)-n]
                list2 = np.array(list2)
                list2 = list2.tolist()
        else:
            list2 = list2
        return list2

    # Finds added pressure at injector due to hydrostatic force, affected by G-forces
    def get_dP4x(self, Fluid, kg, rho, g):
        # Find prop volume from mass and density
        V = kg / rho

        if Fluid == "Oxygen":
            # height of prop in top dome, then adds height of cyl, bottom dome, and lox tank for global x axis coords
            if (V - self.V_LOX_Runline) > (self.V_LOX_cyl + self.V_LOX_dome):
                # WIP
                x = 0
                print('ERROR: WIP PROP LOADED IN TOP LOX DOME')

            # height of prop in tank cylinder, then converted to global x axis coordinates
            elif (V - self.V_LOX_Runline) > self.V_LOX_dome:
                x = (V - self.V_LOX_dome) / (self.pi * self.r_LOX_Tank ** 2) + self.x0_LOX_tank + \
                    self.r_LOX_Tank

            # Height of prop in bottom dome, then convert to global x axis
            elif (V - self.V_LOX_Runline) > 0:
                x = (V / (2 / 3 * self.pi)) ** (1 / 3) + self.x0_LOX_tank

            # Tank empty, we didn't fill, good job
            else:
                print('ERROR: NO PROP IN LOX TANK')
                x = self.x0_LOX_tank

        # RP1 Time now, woohooo
        elif Fluid == "PR::kerosene":
            # all the following find initial height of prop in RP1 Tank
            # height of prop in top dome, then adds height of cyl, bottom dome, and rp1 tank for global x axis coords
            if (V - self.V_RP1_Runline) > (self.V_RP1_cyl + self.V_RP1_dome):
                # WIP
                x = 0
                print('ERROR: WIP PROP LOADED IN TOP RP1 DOME')

            # height of prop in tank cylinder, then converted to global x axis coordinates
            elif (V - self.V_RP1_Runline) > self.V_RP1_dome:
                x = (V - self.V_RP1_dome) / (self.pi * self.r_RP1_Tank ** 2) + self.x0_RP1_tank + \
                    self.r_RP1_Tank

            # Height of prop in bottom dome, then convert to global x axis
            elif (V - self.V_RP1_Runline) > 0:
                x = (V / (2 / 3 * self.pi)) ** (1 / 3) + self.x0_RP1_tank

            # Tank empty, we didn't fill, good job
            else:
                print('ERROR: NO PROP IN RP1 TANK')
                x = self.x0_RP1_tank
        else:
            print('ERROR INVALID FLUID IN get_dP4x')
            x = 0

        # Extra pressure due to hydrostatic effects, with current g accleration
        dP_pa = rho * abs(g) * x
        dP = dP_pa * self.pa2psi
        return dP

    # gets vehicle acceleration from net forces
    def get_a4F(self, thrust_lbf, kg_rkt, drag_N):
        F_thrust = thrust_lbf * self.lbf2N
        F_gravity = kg_rkt * self.g0

        F_net = F_thrust - F_gravity - drag_N
        a = F_net / kg_rkt
        return a

    # finds change in velocity and height for net rocket acceleration and time elapsed
    def get_z4a(self, a, vel, z, dT):
        velnext = vel + a * dT
        znext = z + vel * dT + 1/2 * a * dT ** 2
        return znext, velnext

    # Find atmospheric density and pressure for current altitude using lookup tables
    def get_atm4z(self, z):
        i = 0
        idx_low = 0
        while self.atm_z_lookup[i] < z:
            idx_low = i
            i = i + 1
        idx_high = idx_low + 1

        z_low = self.atm_z_lookup[idx_low]
        z_high = self.atm_z_lookup[idx_high]

        rho_low = self.atm_rho_lookup[idx_low]
        rho_high = self.atm_rho_lookup[idx_high]

        P_low = self.atm_P_lookup[idx_low]
        P_high = self.atm_P_lookup[idx_high]

        per = (z - z_low) / (z_high - z_low)

        rho = (rho_high - rho_low) * per + rho_low
        P = (P_high - P_low) * per + P_low
        P = P * self.pa2psi
        return P, rho

    # gets force due to drag in N from atmospheric density and vehicle velocity
    def get_Fdrag(self, rho_atm, vel, Cd, A):
        if vel >= 0:
            Fdrag = 1/2 * rho_atm * vel**2 * Cd * A

        # makes drag negative to counteract acceleration due to gravity on descent
        else:
            Fdrag = -1 * 1 / 2 * rho_atm * vel ** 2 * Cd * A

        return Fdrag

    # T- Fill and Prepress  --------------------------------------------------------------------------------------------

    # Find Tank and ullage volumes in m^3 and initial prop heights (from engine cap height) in m

    def set_init_prop(self):
        # Set ambient atmospheric pressure and density
        self.P0_atm = self.atm_P_lookup[0] * self.pa2psi
        self.rho0_atm = self.atm_rho_lookup[0]

        self.rho0_LOX = self.get_rho(self.T0_LOX, self.MEOP_LOX, "Oxygen")
        self.V0_LOX = self.kg0_LOX / self.rho0_LOX
        self.V0_LOX_ull = self.V_LOX_Tank - self.V0_LOX
        self.V_LOX_dome = 4 / 6 * self.r_LOX_Tank ** 3
        self.V_LOX_cyl = self.V_LOX_Tank - 2 * self.V_LOX_dome

        # all the following find initial height of prop in LOX Tank

        # Tank Overfilled -> bad
        if (self.V0_LOX-self.V_LOX_Runline) > self.V_LOX_Tank:
            print('ERROR: LOX TANK OVERFILLED')
            self.x0_LOX = 0

        # height of prop in top dome, then adds height of cyl, bottom dome, and lox tank for global x axis coords
        elif (self.V0_LOX-self.V_LOX_Runline) > (self.V_LOX_cyl + self.V_LOX_dome):
            # WIP
            self.x0_LOX = 0
            print('ERROR: WIP PROP LOADED IN TOP LOX DOME')

        # height of prop in tank cylinder, then converted to global x axis coordinates
        elif (self.V0_LOX-self.V_LOX_Runline) > self.V_LOX_dome:
            self.x0_LOX = (self.V0_LOX - self.V_LOX_dome) / (self.pi * self.r_LOX_Tank ** 2) + self.x0_LOX_tank + \
                self.r_LOX_Tank

        # Height of prop in bottom dome, then convert to global x axis
        elif (self.V0_LOX-self.V_LOX_Runline) > 0:
            self.x0_LOX = (self.V0_LOX / (2/3 * self.pi)
                           ) ** (1/3) + self.x0_LOX_tank

        # Tank empty, we didn't fill, good job
        else:
            print('ERROR: NO PROP IN LOX TANK')
            self.x0_LOX = self.x0_LOX_tank

        # RP1 Time now, woohooo

        CP.add_fluids_as_JSON("PR", json.dumps([self.fluid_data_RP1]))
        self.rho0_RP1 = self.get_rho(
            self.T0_RP1, self.MEOP_RP1, "PR::kerosene")
        self.V0_RP1 = self.kg0_RP1 / self.rho0_RP1
        self.V0_RP1_ull = self.V_RP1_Tank - self.V0_RP1
        self.V_RP1_dome = 4 / 6 * self.r_RP1_Tank ** 3
        self.V_RP1_cyl = self.V_RP1_Tank - 2 * self.V_RP1_dome

        # all the following find initial height of prop in RP1 Tank

        # Tank Overfilled -> bad
        if (self.V0_RP1-self.V_RP1_Runline) > self.V_RP1_Tank:
            print('ERROR: RP1 TANK OVERFILLED')
            self.x0_RP1 = 0

        # height of prop in top dome, then adds height of cyl, bottom dome, and rp1 tank for global x axis coords
        elif (self.V0_RP1-self.V_RP1_Runline) > (self.V_RP1_cyl + self.V_RP1_dome):
            # WIP
            self.x0_LOX = 0
            print('ERROR: WIP PROP LOADED IN TOP RP1 DOME')

        # height of prop in tank cylinder, then converted to global x axis coordinates
        elif (self.V0_RP1-self.V_RP1_Runline) > self.V_RP1_dome:
            self.x0_RP1 = (self.V0_RP1 - self.V_RP1_dome) / (self.pi * self.r_RP1_Tank ** 2) + self.x0_RP1_tank + \
                self.r_RP1_Tank

        # Height of prop in bottom dome, then convert to global x axis
        elif (self.V0_RP1-self.V_RP1_Runline) > 0:
            self.x0_RP1 = (self.V0_RP1 / (2/3 * self.pi)
                           ) ** (1/3) + self.x0_RP1_tank

        # Tank empty, we didn't fill, good job
        else:
            print('ERROR: NO PROP IN RP1 TANK')
            self.x0_RP1 = self.x0_RP1_tank

        # Initializes mass of helium in tank
        self.kg0_GHE = self.get_kg4dP(
            self.V_GHE_Tank, self.MEOP_GHE, self.T0_GHE, 1, self.fluidPress)

        # Set initial wet mass of rocket
        self.kg0_rkt = self.kg0_GHE + self.kg0_RP1 + self.kg0_LOX + self.kg_drymass

        # Set initial hydrostatic pressures
        self.P0_LOX_hydro = self.get_dP4x(
            'Oxygen', self.kg0_LOX, self.rho0_LOX, self.g0)
        self.P0_RP1_hydro = self.get_dP4x(
            'PR::kerosene', self.kg0_RP1, self.rho0_RP1, self.g0)

        # set initial acceleration of rocket (-9.81m/s)
        self.a0 = self.get_a4F(0, self.kg0_rkt, 0)

    # finds lumped CdAs for runlines using ideal conditions
    def set_propCdA(self):
        mdot_LOX = self.mdot_LOX_const
        mdot_RP1 = self.mdot_RP1_const

        P_LOX_tank = self.MEOP_LOX + self.P0_LOX_hydro
        P_RP1_tank = self.MEOP_RP1 + self.P0_RP1_hydro

        P_mid_LOX = (P_LOX_tank + self.MEOP_Pc)/2
        P_mid_RP1 = (P_RP1_tank + self.MEOP_Pc)/2

        rho_LOX = self.get_rho(self.T0_LOX, P_mid_LOX, "Oxygen")
        rho_RP1 = self.get_rho(self.T0_RP1, P_mid_RP1, "PR::kerosene")

        dP_LOX = (P_LOX_tank - self.MEOP_Pc) * self.psi2pa
        dP_RP1 = (P_RP1_tank - self.MEOP_Pc) * self.psi2pa

        self.CdA_LOX = mdot_LOX / np.sqrt(2 * rho_LOX * dP_LOX)  # CdA in m^2
        self.CdA_RP1 = mdot_RP1 / np.sqrt(2 * rho_RP1 * dP_RP1)  # CdA in m^2

    def prepress(self):
        """ Trajectory Markers
        print('- Prepress')
        print('  * Begin Prepress')
        """

        self.P_GHE[0] = self.MEOP_GHE
        self.T_GHE[0] = self.T0_GHE
        self.P_LOX[0] = self.P0_atm
        self.P_RP1[0] = self.P0_atm
        self.kg_GHE[0] = self.kg0_GHE
        # WIP needs to be at kg for voume at Pamb
        self.kg_LOX_ull = [0] * len(self.time)
        self.kg_RP1_ull = [0] * len(self.time)  # WIP ^^

        i = 0
        while (self.P_RP1[i] < self.Pmax_RP1) | (self.P_LOX[i] < self.Pmax_LOX):
            self.time.append(0)
            self.P_GHE.append(0)
            self.T_GHE.append(0)
            self.P_RP1.append(0)
            self.P_LOX.append(0)
            self.kg_GHE.append(0)
            self.kg_LOX_ull.append(0)
            self.kg_RP1_ull.append(0)

            # Looks at current tank pressures and bang bang flags to decide the next iteration bang bang states
            self.set_pressConfig(self.P_LOX[i], self.P_RP1[i], 0)

            # Uses bang bang config with current pressures to find # of kg moving to LOX, RP1 tank and RCS, leaving GHE
            dkg_LOX, dkg_RP1, dkg_RCS, dkg_GHE = self.get_pressdKg(
                self.P_LOX[i], self.P_RP1[i], self.P0_atm, self.P_GHE[i], self.T_GHE[i])

            # Find ullage helium kg count
            kg_LOX_ull_next = self.kg_LOX_ull[i] + dkg_LOX
            kg_RP1_ull_next = self.kg_RP1_ull[i] + dkg_RP1

            # Output GHE kgs
            kg_GHE_next = self.kg_GHE[i] + dkg_GHE

            # Uses dkg info to find new tank equalization pressures and temperatures
            P_LOX_next, P_RP1_next, P_GHE_next, T_GHE_next = self.get_pressEq(
                self.V0_LOX_ull, kg_LOX_ull_next, self.V0_RP1_ull, kg_RP1_ull_next, self.T_GHE[i], kg_GHE_next)

            time_next = self.time[i] + self.dT

            # Iterates index and defines the next index properties
            i = i + 1
            self.P_GHE[i] = P_GHE_next
            self.T_GHE[i] = T_GHE_next
            self.P_LOX[i] = P_LOX_next
            self.P_RP1[i] = P_RP1_next
            self.kg_GHE[i] = kg_GHE_next
            self.time[i] = time_next
            self.kg_LOX_ull[i] = kg_LOX_ull_next
            self.kg_RP1_ull[i] = kg_RP1_ull_next

            """ Trajectory Markers 
        print('  * Prepress Complete')
        print('')

        print('- At-Pressure Hold')
        print('  * Begin Hold')
            """

        hold = self.prepress_hold
        imax = int(hold / self.dT)

        P_GHE = self.P_GHE[i]
        T_GHE = self.T_GHE[i]
        P_LOX = self.P_LOX[i]
        P_RP1 = self.P_RP1[i]
        kg_GHE = self.kg_GHE[i]
        time = self.time[i]
        kg_LOX_ull = self.kg_LOX_ull[i]
        kg_RP1_ull = self.kg_RP1_ull[i]

        i = 1
        for i in range(i, imax):
            self.P_GHE.append(P_GHE)
            self.T_GHE.append(T_GHE)
            self.P_LOX.append(P_LOX)
            self.P_RP1.append(P_RP1)
            self.kg_GHE.append(kg_GHE)
            self.time.append(time + i * self.dT)
            self.kg_LOX_ull.append(kg_LOX_ull)
            self.kg_RP1_ull.append(kg_RP1_ull)

            """ Trajectory Markers 
        print('  * Hold Complete')
        print('')
            """

    # Burn Baby Burn (static fire or launch)  --------------------------------------------------------------------------

    # def staticfire(self):
    #     # Initialize new variables that were not needed during prepress
    #     self.kg_RP1 = [self.kg0_RP1] * len(self.time)
    #     self.kg_LOX = [self.kg0_LOX] * len(self.time)
    #     self.V_LOX_ull = [self.V0_LOX_ull] * len(self.time)
    #     self.V_RP1_ull = [self.V0_RP1_ull] * len(self.time)
    #     self.V_LOX = [self.V0_LOX] * len(self.time)
    #     self.V_RP1 = [self.V0_RP1] * len(self.time)
    #
    #     self.Pc = [0] * len(self.time)
    #     self.mdot_LOX = [0] * len(self.time)
    #     self.mdot_RP1 = [0] * len(self.time)
    #     self.OtF = [0] * len(self.time)
    #     self.gamma = [0] * len(self.time)
    #     self.gasConst = [0] * len(self.time)
    #     self.T_flame = [0] * len(self.time)
    #     self.thrust = [0] * len(self.time)
    #     self.Isp = [0] * len(self.time)
    #
    #     self.kg_rkt = [self.kg0_rkt] * len(self.time)
    #
    #     self.rho_LOX = [self.rho0_LOX] * len(self.time)
    #     self.rho_RP1 = [self.rho0_RP1] * len(self.time)
    #
    #     self.P_LOX_hydro = [self.P0_LOX_hydro] * len(self.time)
    #     self.P_RP1_hydro = [self.P0_RP1_hydro] * len(self.time)
    #
    #     # I index starting at end of prepress
    #     i = len(self.time) - 1
    #     # Initializing First combustion parameters
    #     self.OtF[i] = self.mdot_LOX_const / self.mdot_RP1_const
    #     self.Pc[i] = self.MEOP_Pc
    #
    #     # Stops when growing ullage volume is equal to COPV volume
    #     print('- Ignition!')
    #     n = 1 - self.n_shutoff # shut off with ##% volume prop in tanks as margin
    #     while (self.V_LOX_ull[i] < self.V_LOX_Tank * n) & (self.V_RP1_ull[i] < self.V_RP1_Tank * n):
    #         self.time.append(0)
    #         self.P_GHE.append(0)
    #         self.T_GHE.append(0)
    #         self.P_RP1.append(0)
    #         self.P_LOX.append(0)
    #         self.kg_GHE.append(0)
    #         self.kg_LOX.append(0)
    #         self.kg_RP1.append(0)
    #         self.V_LOX_ull.append(0)
    #         self.V_RP1_ull.append(0)
    #         self.kg_LOX_ull.append(0)
    #         self.kg_RP1_ull.append(0)
    #         self.V_LOX.append(0)
    #         self.V_RP1.append(0)
    #         self.kg_rkt.append(0)
    #         # Non steady prop mdots
    #         self.Pc.append(0)
    #         self.OtF.append(0)
    #         self.gamma.append(0)
    #         self.gasConst.append(0)
    #         self.T_flame.append(0)
    #         self.thrust.append(0)
    #         self.Isp.append(0)
    #
    #         self.rho_LOX.append(0)
    #         self.rho_RP1.append(0)
    #
    #         self.P_LOX_hydro.append(0)
    #         self.P_RP1_hydro.append(0)
    #
    #         # Looks at current tank pressures and bang bang flags to decide the next iteration bang bang states
    #         self.set_pressConfig(self.P_LOX[i], self.P_RP1[i], 0)
    #
    #         # Uses bang bang config with current pressures to find # of kg moving to LOX, RP1 tank and RCS, leaving GHE
    #         dkg_LOX, dkg_RP1, dkg_RCS, dkg_GHE = self.get_pressdKg(self.P_LOX[i], self.P_RP1[i], self.P0_atm,
    #                                                                self.P_GHE[i], self.T_GHE[i])
    #
    #         # Find ullage helium kg count
    #         kg_LOX_ull_next = self.kg_LOX_ull[i] + dkg_LOX
    #         kg_RP1_ull_next = self.kg_RP1_ull[i] + dkg_RP1
    #
    #         # Output GHE kgs
    #         kg_GHE_next = self.kg_GHE[i] + dkg_GHE
    #
    #         # Uses dkg info to find new tank equalization pressures and temperatures
    #         P_LOX_next, P_RP1_next, P_GHE_next, T_GHE_next = self.get_pressEq(self.V_LOX_ull[i], kg_LOX_ull_next,
    #                                                                           self.V_RP1_ull[i], kg_RP1_ull_next,
    #                                                                           self.T_GHE[i], kg_GHE_next)
    #
    #         P_LOX_tot = self.P_LOX[i] + self.P_LOX_hydro[i]
    #         P_RP1_tot = self.P_RP1[i] + self.P_RP1_hydro[i]
    #
    #         # Find prop mdots for given chamber pressure
    #         self.mdot_LOX[i] = self.get_propMdot(self.Pc[i], P_LOX_tot, "Oxygen")
    #         self.mdot_RP1[i] = self.get_propMdot(self.Pc[i], P_RP1_tot, "PR::kerosene")
    #
    #         # Moving Average for incoming mass flow, acts as dampening on mdot oscillations wooo (moving mean, movmean)
    #         mdot_LOX_last = self.get_movmean(self.mdot_LOX, self.n_movmean)
    #         self.mdot_LOX[i] = mdot_LOX_last[len(mdot_LOX_last) - 1]
    #
    #         mdot_RP1_last = self.get_movmean(self.mdot_RP1, self.n_movmean)
    #         self.mdot_RP1[i] = mdot_RP1_last[len(mdot_RP1_last) - 1]
    #
    #         # appends had to be moved below the realtime moving average for some reason, i dont know why, I do not claim to understand the python gods
    #         self.mdot_LOX.append(0)
    #         self.mdot_RP1.append(0)
    #
    #         # Find delta in ullage volume
    #         V_LOX_ull_next, V_RP1_ull_next = self.get_dVullage(self.mdot_LOX[i], self.P_LOX[i], self.T0_LOX, self.V_LOX_ull[i], self.mdot_RP1[i], self.P_RP1[i], self.T0_RP1, self.V_RP1_ull[i])
    #
    #         time_next = self.time[i] + self.dT
    #
    #         # iterating i to next index, start solving for inlet conditions of next time index
    #         i = i + 1
    #         self.V_LOX_ull[i] = V_LOX_ull_next
    #         self.V_RP1_ull[i] = V_RP1_ull_next
    #         self.P_GHE[i] = P_GHE_next
    #         self.T_GHE[i] = T_GHE_next
    #         self.P_LOX[i] = P_LOX_next
    #         self.P_RP1[i] = P_RP1_next
    #         self.kg_GHE[i] = kg_GHE_next
    #         self.time[i] = time_next
    #         self.kg_LOX_ull[i] = kg_LOX_ull_next
    #         self.kg_RP1_ull[i] = kg_RP1_ull_next
    #
    #         self.rho_LOX[i] = self.get_rho(self.T0_LOX, self.P_LOX[i], 'Oxygen')
    #         self.rho_RP1[i] = self.get_rho(self.T0_RP1, self.P_RP1[i], 'PR::kerosene')
    #
    #         # COMBUSTION TIME WOOOOO
    #         self.OtF[i] = self.mdot_LOX[i-1] / self.mdot_RP1[i-1]
    #
    #         # Get gamma, R, and flame T from OtF and Pc
    #         self.gamma[i], self.gasConst[i], self.T_flame[i] = self.get_CEA(self.OtF[i], self.Pc[i-1])
    #
    #         # Get new Pc from R flame T and gamma throat area and prop mdots
    #         self.Pc[i] = self.get_Pc4CEA(self.gamma[i], self.gasConst[i], self.T_flame[i], (self.mdot_LOX[i-1]+self.mdot_RP1[i-1]))
    #
    #         mdot_tot = self.mdot_LOX[i-1]+self.mdot_RP1[i-1]
    #
    #         # Get thrust from CEA values
    #         self.thrust[i], self.Isp[i] = self.get_thrust4CEA(self.gamma[i], self.gasConst[i], self.T_flame[i], self.P0_atm, self.Pc[i], mdot_tot)
    #
    #
    #         # Find new prop tank pressures for expanded ullage volume
    #         self.P_LOX[i] = self.get_P4V(V_LOX_ull_next, kg_LOX_ull_next, self.T_GHE[i], 2, self.fluidPress)
    #         self.P_RP1[i] = self.get_P4V(V_RP1_ull_next, kg_RP1_ull_next, self.T_GHE[i], 1, self.fluidPress)
    #
    #         # Kg change in prop due to engine yeeting it out the rocket
    #         self.kg_RP1[i] = self.kg_RP1[i - 1] - self.mdot_RP1[i-1] * self.dT
    #         self.kg_LOX[i] = self.kg_LOX[i - 1] - self.mdot_LOX[i-1] * self.dT
    #
    #         # Vol change in propellants
    #         self.V_LOX[i] = self.V_LOX[i - 1] - (V_LOX_ull_next - self.V_LOX_ull[i-1])
    #         self.V_RP1[i] = self.V_RP1[i - 1] - (V_RP1_ull_next - self.V_RP1_ull[i-1])
    #
    #         # Change in rocket total mass
    #         self.kg_rkt[i] = self.kg_rkt[i-1] - (self.mdot_RP1[i-1] + self.mdot_LOX[i-1]) * self.dT
    #
    #         # Change in rocket propellant hydrostatic pressure
    #         self.P_LOX_hydro[i] = self.get_dP4x('Oxygen', self.kg_LOX[i], self.rho_LOX[i], self.g0)
    #         self.P_RP1_hydro[i] = self.get_dP4x('PR::kerosene', self.kg_RP1[i], self.rho_RP1[i], self.g0)
    #
    #
    #     print('- MECO, Prop Tanks Depleted')
    #     # Mass Fill Percentage
    #     self.kg_GHE_per = [(x / self.kg0_GHE) * 100 for x in self.kg_GHE]
    #     self.kg_RP1_per = [(x / self.kg0_RP1) * 100 for x in self.kg_RP1]
    #     self.kg_LOX_per = [(x / self.kg0_LOX) * 100 for x in self.kg_LOX]
    #
    #     # Volume Fill Percentage
    #     self.V_RP1_per = [(x / self.V_RP1_Tank) * 100 for x in self.V_RP1]
    #     self.V_LOX_per = [(x / self.V_LOX_Tank) * 100 for x in self.V_LOX]

    def launch(self):
        # Initialize new variables that were not needed during prepress
        self.kg_RP1 = [self.kg0_RP1] * len(self.time)
        self.kg_LOX = [self.kg0_LOX] * len(self.time)
        self.V_LOX_ull = [self.V0_LOX_ull] * len(self.time)
        self.V_RP1_ull = [self.V0_RP1_ull] * len(self.time)
        self.V_LOX = [self.V0_LOX] * len(self.time)
        self.V_RP1 = [self.V0_RP1] * len(self.time)

        self.Pc = [0] * len(self.time)
        self.mdot_LOX = [0] * len(self.time)
        self.mdot_RP1 = [0] * len(self.time)
        self.OtF = [0] * len(self.time)
        self.gamma = [0] * len(self.time)
        self.gasConst = [0] * len(self.time)
        self.T_flame = [0] * len(self.time)
        self.thrust = [0] * len(self.time)
        self.Isp = [0] * len(self.time)

        self.kg_rkt = [self.kg0_rkt] * len(self.time)

        self.rho_LOX = [self.rho0_LOX] * len(self.time)
        self.rho_RP1 = [self.rho0_RP1] * len(self.time)

        self.P_LOX_hydro = [self.P0_LOX_hydro] * len(self.time)
        self.P_RP1_hydro = [self.P0_RP1_hydro] * len(self.time)

        self.a = [self.a0] * len(self.time)
        self.vel = [0] * len(self.time)
        self.z = [0] * len(self.time)

        self.P_atm = [self.P0_atm] * len(self.time)
        self.rho_atm = [self.rho0_atm] * len(self.time)

        self.Fdrag = [0] * len(self.time)

        self.impulse = [0] * len(self.time)

        # I index starting at end of prepress
        i = len(self.time) - 1
        # Initializing First combustion parameters
        self.T_0 = self.time[i]
        self.OtF[i] = self.mdot_LOX_const / self.mdot_RP1_const
        self.Pc[i] = self.MEOP_Pc

        # Stops when growing ullage volume is equal to COPV volume

        """ Trajectory Markers 
        print('- Ignition')
        print('')
        """
        n = 1 - self.n_shutoff  # shut off with ##% volume prop in tanks as margin
        while (self.V_LOX_ull[i] < self.V_LOX_Tank * n) & (self.V_RP1_ull[i] < self.V_RP1_Tank * n) & (self.impulse[i] < self.maxImpulse):
            self.time.append(0)
            self.P_GHE.append(0)
            self.T_GHE.append(0)
            self.P_RP1.append(0)
            self.P_LOX.append(0)
            self.kg_GHE.append(0)
            self.kg_LOX.append(0)
            self.kg_RP1.append(0)
            self.V_LOX_ull.append(0)
            self.V_RP1_ull.append(0)
            self.kg_LOX_ull.append(0)
            self.kg_RP1_ull.append(0)
            self.V_LOX.append(0)
            self.V_RP1.append(0)
            self.kg_rkt.append(0)
            # Non steady prop mdots
            self.Pc.append(0)
            self.OtF.append(0)
            self.gamma.append(0)
            self.gasConst.append(0)
            self.T_flame.append(0)
            self.thrust.append(0)
            self.Isp.append(0)

            self.rho_LOX.append(0)
            self.rho_RP1.append(0)

            self.P_LOX_hydro.append(0)
            self.P_RP1_hydro.append(0)

            self.a.append(0)
            self.vel.append(0)
            self.z.append(0)

            self.P_atm.append(0)
            self.rho_atm.append(0)

            self.Fdrag.append(0)

            self.impulse.append(0)

            # Looks at current tank pressures and bang bang flags to decide the next iteration bang bang states
            self.set_pressConfig(self.P_LOX[i], self.P_RP1[i], 0)

            # Uses bang bang config with current pressures to find # of kg moving to LOX, RP1 tank and RCS, leaving GHE
            dkg_LOX, dkg_RP1, dkg_RCS, dkg_GHE = self.get_pressdKg(self.P_LOX[i], self.P_RP1[i], self.P_atm[i],
                                                                   self.P_GHE[i], self.T_GHE[i])

            # Find ullage helium kg count
            kg_LOX_ull_next = self.kg_LOX_ull[i] + dkg_LOX
            kg_RP1_ull_next = self.kg_RP1_ull[i] + dkg_RP1

            # Output GHE kgs
            kg_GHE_next = self.kg_GHE[i] + dkg_GHE

            # Uses dkg info to find new tank equalization pressures and temperatures
            P_LOX_next, P_RP1_next, P_GHE_next, T_GHE_next = self.get_pressEq(self.V_LOX_ull[i], kg_LOX_ull_next,
                                                                              self.V_RP1_ull[i], kg_RP1_ull_next,
                                                                              self.T_GHE[i], kg_GHE_next)

            P_LOX_tot = self.P_LOX[i] + self.P_LOX_hydro[i]
            P_RP1_tot = self.P_RP1[i] + self.P_RP1_hydro[i]

            # Find prop mdots for given chamber pressure
            self.mdot_LOX[i] = self.get_propMdot(
                self.Pc[i], P_LOX_tot, "Oxygen")
            self.mdot_RP1[i] = self.get_propMdot(
                self.Pc[i], P_RP1_tot, "PR::kerosene")

            # Moving Average for incoming mass flow, acts as dampening on mdot oscillations wooo (moving mean, movmean)
            mdot_LOX_last = self.get_movmean(self.mdot_LOX, self.n_movmean)
            self.mdot_LOX[i] = mdot_LOX_last[len(mdot_LOX_last) - 1]

            mdot_RP1_last = self.get_movmean(self.mdot_RP1, self.n_movmean)
            self.mdot_RP1[i] = mdot_RP1_last[len(mdot_RP1_last) - 1]

            # appends had to be moved below the realtime moving average for some reason, i dont know why, I do not claim to understand the python gods
            self.mdot_LOX.append(0)
            self.mdot_RP1.append(0)

            # Find delta in ullage volume
            V_LOX_ull_next, V_RP1_ull_next = self.get_dVullage(
                self.mdot_LOX[i], self.P_LOX[i], self.T0_LOX, self.V_LOX_ull[i], self.mdot_RP1[i], self.P_RP1[i], self.T0_RP1, self.V_RP1_ull[i])

            time_next = self.time[i] + self.dT

            # iterating i to next index, start solving for inlet conditions of next time index
            i = i + 1

            self.V_LOX_ull[i] = V_LOX_ull_next
            self.V_RP1_ull[i] = V_RP1_ull_next
            self.P_GHE[i] = P_GHE_next
            self.T_GHE[i] = T_GHE_next
            self.P_LOX[i] = P_LOX_next
            self.P_RP1[i] = P_RP1_next
            self.kg_GHE[i] = kg_GHE_next
            self.time[i] = time_next
            self.kg_LOX_ull[i] = kg_LOX_ull_next
            self.kg_RP1_ull[i] = kg_RP1_ull_next

            self.rho_LOX[i] = self.get_rho(
                self.T0_LOX, self.P_LOX[i], 'Oxygen')
            self.rho_RP1[i] = self.get_rho(
                self.T0_RP1, self.P_RP1[i], 'PR::kerosene')

            # COMBUSTION TIME WOOOOO
            self.OtF[i] = self.mdot_LOX[i-1] / self.mdot_RP1[i-1]

            # Get gamma, R, and flame T from OtF and Pc
            self.gamma[i], self.gasConst[i], self.T_flame[i] = self.get_CEA(
                self.OtF[i], self.Pc[i-1])

            # Get new Pc from R flame T and gamma throat area and prop mdots
            self.Pc[i] = self.get_Pc4CEA(
                self.gamma[i], self.gasConst[i], self.T_flame[i], (self.mdot_LOX[i-1]+self.mdot_RP1[i-1]))

            mdot_tot = (self.mdot_LOX[i-1]+self.mdot_RP1[i-1])

            # Get thrust from CEA values
            self.thrust[i], self.Isp[i] = self.get_thrust4CEA(
                self.gamma[i], self.gasConst[i], self.T_flame[i], self.P_atm[i-1], self.Pc[i], mdot_tot)

            # Calculates total impulse
            self.impulse[i] = self.impulse[i-1] + self.thrust[i] * self.dT

            # Find new prop tank pressures for expanded ullage volume
            self.P_LOX[i] = self.get_P4V(
                V_LOX_ull_next, kg_LOX_ull_next, self.T_GHE[i], 2, self.fluidPress)
            self.P_RP1[i] = self.get_P4V(
                V_RP1_ull_next, kg_RP1_ull_next, self.T_GHE[i], 1, self.fluidPress)

            # Kg change in prop due to engine yeeting it out the rocket
            self.kg_RP1[i] = self.kg_RP1[i - 1] - self.mdot_RP1[i-1] * self.dT
            self.kg_LOX[i] = self.kg_LOX[i - 1] - self.mdot_LOX[i-1] * self.dT

            # Vol change in propellants
            self.V_LOX[i] = self.V_LOX[i - 1] - \
                (V_LOX_ull_next - self.V_LOX_ull[i-1])
            self.V_RP1[i] = self.V_RP1[i - 1] - \
                (V_RP1_ull_next - self.V_RP1_ull[i-1])

            # Change in rocket total mass
            self.kg_rkt[i] = self.kg_rkt[i-1] - \
                (self.mdot_RP1[i-1] + self.mdot_LOX[i-1]) * self.dT

            # Change in rocket propellant hydrostatic pressure
            self.P_LOX_hydro[i] = self.get_dP4x(
                'Oxygen', self.kg_LOX[i], self.rho_LOX[i], self.a[i-1])
            self.P_RP1_hydro[i] = self.get_dP4x(
                'PR::kerosene', self.kg_RP1[i], self.rho_RP1[i], self.a[i-1])

            # find drag force from previous iteration
            self.Fdrag[i] = self.get_Fdrag(
                self.rho_atm[i-1], self.vel[i-1], 0.32022, (16/2 * self.in2m)**2 * self.pi)

            # find change in acceleration, velocity, and height for given thrust
            self.a[i] = self.get_a4F(
                self.thrust[i], self.kg_rkt[i], self.Fdrag[i])
            self.z[i], self.vel[i] = self.get_z4a(
                self.a[i], self.vel[i-1], self.z[i-1], self.dT)

            self.P_atm[i], self.rho_atm[i] = self.get_atm4z(self.z[i])

            """ Trajectory Markers
        print('- MECO, Prop Tanks Depleted')
        print('  * Impulse of', self.impulse[i], 'lbf-s reached')
        if self.impulse[i] > self.maxImpulse:
            print('  * SHUTDOWN: Max Impulse Bound')
            print('  * Residual LOX: ',
                  self.kg_LOX[i] - (self.V_LOX_Runline) * self.rho_LOX[i], 'kgs')
            print('  * Residual RP1: ',
                  self.kg_RP1[i] - (self.V_RP1_Runline) * self.rho_RP1[i], 'kgs')
            """

        self.MECO = self.time[i]
        burntime = '%.1f' % (self.MECO - self.T_0)

        """ Trajectory Markers 
        print('  *', burntime, 's Burn')
        print('')

        print('- Coast')
        """

        # Post MECO coast
        while self.z[i] > 0:
            self.time.append(0)
            self.P_GHE.append(0)
            self.T_GHE.append(0)
            self.P_RP1.append(0)
            self.P_LOX.append(0)
            self.kg_GHE.append(0)
            self.kg_LOX.append(0)
            self.kg_RP1.append(0)
            self.V_LOX_ull.append(0)
            self.V_RP1_ull.append(0)
            self.kg_LOX_ull.append(0)
            self.kg_RP1_ull.append(0)
            self.V_LOX.append(0)
            self.V_RP1.append(0)
            self.kg_rkt.append(0)
            # Non steady prop mdots
            self.Pc.append(0)
            self.OtF.append(0)
            self.gamma.append(0)
            self.gasConst.append(0)
            self.T_flame.append(0)
            self.thrust.append(0)
            self.Isp.append(0)

            self.rho_LOX.append(0)
            self.rho_RP1.append(0)

            self.mdot_LOX.append(0)
            self.mdot_RP1.append(0)

            self.P_LOX_hydro.append(0)
            self.P_RP1_hydro.append(0)

            self.a.append(0)
            self.vel.append(0)
            self.z.append(0)

            self.P_atm.append(0)
            self.rho_atm.append(0)

            self.Fdrag.append(0)

            time_next = self.time[i] + self.dT

            i = i + 1
            self.V_LOX_ull[i] = self.V_LOX_ull[i-1]
            self.V_RP1_ull[i] = self.V_LOX_ull[i-1]
            self.P_GHE[i] = self.P_GHE[i-1]
            self.T_GHE[i] = self.T_GHE[i-1]
            self.P_LOX[i] = self.P_LOX[i-1]
            self.P_RP1[i] = self.P_RP1[i-1]
            self.kg_GHE[i] = self.kg_GHE[i-1]
            self.time[i] = time_next
            self.kg_LOX_ull[i] = self.kg_LOX_ull[i-1]
            self.kg_RP1_ull[i] = self.kg_LOX_ull[i-1]

            self.rho_LOX[i] = self.get_rho(
                self.T0_LOX, self.P_LOX[i], 'Oxygen')
            self.rho_RP1[i] = self.get_rho(
                self.T0_RP1, self.P_RP1[i], 'PR::kerosene')

            # find drag force from previous iteration
            self.Fdrag[i] = self.get_Fdrag(self.rho_atm[i - 1], self.vel[i - 1], 0.32022,
                                           (16 / 2 * self.in2m) ** 2 * self.pi)

            # COMBUSTION TIME WOOOOO
            self.mdot_LOX[i] = 0
            self.mdot_RP1[i] = 0

            self.OtF[i] = 0

            # Get gamma, R, and flame T from OtF and Pc
            self.gamma[i], self.gasConst[i], self.T_flame[i] = 0, 0, 0

            # Get new Pc from R flame T and gamma throat area and prop mdots
            self.Pc[i] = 0

            # Get thrust from CEA values
            self.thrust[i], self.Isp[i] = 0, 0

            # Find new prop tank pressures for expanded ullage volume
            self.P_LOX[i] = self.P_LOX[i-1]
            self.P_RP1[i] = self.P_RP1[i-1]

            # Kg change in prop due to engine yeeting it out the rocket
            self.kg_RP1[i] = self.kg_RP1[i-1]
            self.kg_LOX[i] = self.kg_LOX[i-1]

            # Vol change in propellants
            self.V_LOX[i] = self.V_LOX[i - 1]
            self.V_RP1[i] = self.V_RP1[i - 1]

            # Change in rocket total mass
            self.kg_rkt[i] = self.kg_rkt[i - 1]

            # Change in rocket propellant hydrostatic pressure
            self.P_LOX_hydro[i] = self.get_dP4x(
                'Oxygen', self.kg_LOX[i], self.rho_LOX[i], self.a[i - 1])
            self.P_RP1_hydro[i] = self.get_dP4x(
                'PR::kerosene', self.kg_RP1[i], self.rho_RP1[i], self.a[i - 1])

            # find change in acceleration, velocity, and height for given thrust
            self.a[i] = self.get_a4F(
                self.thrust[i], self.kg_rkt[i], self.Fdrag[i])
            self.z[i], self.vel[i] = self.get_z4a(
                self.a[i], self.vel[i - 1], self.z[i - 1], self.dT)

            # Finds new atmospheric conditions for given z
            self.P_atm[i], self.rho_atm[i] = self.get_atm4z(self.z[i])

            """ Trajectory Markers 
        print('  * Apogee of', max(self.z)/1000, 'km Reached')
        print('')
        print('- Ballistic Touchdown')
            """

        # Mass Fill Percentage
        self.kg_GHE_per = [(x / self.kg0_GHE) * 100 for x in self.kg_GHE]
        self.kg_RP1_per = [(x / self.kg0_RP1) * 100 for x in self.kg_RP1]
        self.kg_LOX_per = [(x / self.kg0_LOX) * 100 for x in self.kg_LOX]

        # Volume Fill Percentage
        self.V_RP1_per = [(x / self.V_RP1_Tank) * 100 for x in self.V_RP1]
        self.V_LOX_per = [(x / self.V_LOX_Tank) * 100 for x in self.V_LOX]

        # Data Transfer to MATLAB
        mdot_LOX_RP1_pers.append([self.time, self.mdot_RP1, self.mdot_LOX])

# ----------------------------------------------------------------------------------------------------------------------

    # Lookin' For Some Plot Stuff Baby This Evening

    # Just normal plots
    def plot(self):
        pl.figure(1)

        pl.subplot(2, 3, 1)
        pl.title('Helium Tank Pressure')
        pl.plot(self.time, self.P_GHE, 'm', label='Helium Tank Pressure')
        pl.xlabel('Time (s)')
        pl.ylabel('Pressure (psia)')
        pl.grid()
        pl.legend()

        pl.subplot(2, 3, 2)
        pl.title('Propellant Tank Pressures')
        Pmax_LOX = [self.Pmax_LOX] * len(self.time)
        Pmin_LOX = [self.Pmin_LOX] * len(self.time)

        Pmax_RP1 = [self.Pmax_RP1] * len(self.time)
        Pmin_RP1 = [self.Pmin_RP1] * len(self.time)

        pl.plot(self.time, Pmax_LOX, ':b')
        pl.plot(self.time, Pmin_LOX, ':b')
        pl.plot(self.time, self.P_LOX, 'b', label='LOX Tank Pressure')
        pl.plot(self.time, self.P_RP1, 'r', label='RP1 Tank Pressure')
        pl.plot(self.time, Pmax_RP1, ':r')
        pl.plot(self.time, Pmin_RP1, ':r')
        pl.xlabel('Time (s)')
        pl.ylabel('Pressure (psia)')
        pl.grid()
        pl.legend()

        # pl.subplot(2,3,3)
        # pl.title('Propellant Tank Ullage Mass')
        # pl.plot(self.time, self.kg_GHE, 'm', label='Helium Press Mass')
        # pl.plot(self.time, self.kg_LOX_ull, 'b', label='LOX Ullage Mass')
        # pl.plot(self.time, self.kg_RP1_ull, 'r', label='RP1 Ullage Mass')
        # pl.xlabel('Time (s)')
        # pl.ylabel('Mass (kg)')
        # pl.grid()
        # pl.legend()

        # pl.subplot(2, 3, 4)
        # pl.title('Helium Tank Temperature')
        # pl.plot(self.time, self.T_GHE, 'm', label='Helium Tank Temperature')
        # pl.xlabel('Time (s)')
        # pl.ylabel('Temperature (K)')
        # pl.grid()
        # pl.legend()

        pl.subplot(2, 3, 3)
        pl.title('Propellant Mass Flow Rate')
        pl.plot(self.time, self.mdot_LOX, 'b', label='mdot LOX (kg/s)')
        pl.plot(self.time, self.mdot_RP1, 'r', label='mdot RP1 (kg/s)')
        pl.xlabel('Time (s)')
        pl.ylabel('Mass Flow Rate (kg/s)')
        pl.grid()
        pl.legend()

        # pl.subplot(2, 3, 4)
        # pl.title('Propellant Tank Volume Percent Fill')
        # pl.plot(self.time, self.V_RP1_per, 'r', label='RP1 Tank Volume % Fill')
        # pl.plot(self.time, self.V_LOX_per, 'b', label='LOX Tank Volume % Fill')
        # pl.xlabel('Time (s)')
        # pl.ylabel('Volume % Fill')
        # pl.grid()
        # pl.legend()

        pl.subplot(2, 3, 4)
        pl.title('Propellant Tank Mass Percent Fill')
        pl.plot(self.time, self.kg_GHE_per, 'm',
                label='Helium Tank Mass Percent Fill')
        pl.plot(self.time, self.kg_RP1_per, 'r',
                label='RP1 Tank Mass Percent Fill')
        pl.plot(self.time, self.kg_LOX_per, 'b',
                label='LOX Tank Mass Percent Fill')
        pl.xlabel('Time (s)')
        pl.ylabel('Mass Percent Fill')
        pl.grid()
        pl.legend()

        # pl.subplot(2, 3, 5)
        # pl.title('Chamber Pressure')
        # pl.plot(self.time, self.Pc, 'orange', label='Chamber Pressure')
        # OtFx100 = [x * 100 for x in self.OtF]
        # Pcnom = [self.MEOP_Pc] * len(self.time)
        # OtFnom = [self.mdot_LOX_const /
        #          self.mdot_RP1_const * 100] * len(self.time)
        # pl.plot(self.time, Pcnom, '--', color='orange')
        # pl.plot(self.time, OtFnom, '--g')
        # pl.plot(self.time, OtFx100, 'g', label='O/F Ratio x 100')
        # pl.xlabel('Time (s)')
        # pl.ylabel('Pressure (psia)')
        # pl.grid()
        # pl.legend()

        # pl.subplot(2, 3, 6)
        # pl.title('Engine Thrust')
        # pl.plot(self.time, self.thrust, 'hotpink', label='Thrust (lbf)')
        # pl.plot(self.time, self.Isp, 'c', label='Specific Impulse (s)')
        # thrustnom = [3500] * len(self.time)
        # pl.plot(self.time, thrustnom, '--', color='hotpink')
        # pl.xlabel('Time (s)')
        # pl.ylabel('Thrust (lbf) / Isp (s)')
        # pl.grid()
        # pl.legend()

        print('- Static Plot Complete')
        pl.show()

        """

        pl.figure(2)

        pl.subplot(2, 3, 1)
        pl.title('Chamber Pressure')
        pl.plot(self.time, self.Pc, 'orange', label='Chamber Pressure')
        OtFx100 = [x * 100 for x in self.OtF]
        Pcnom = [self.MEOP_Pc] * len(self.time)
        OtFnom = [self.mdot_LOX_const /
                  self.mdot_RP1_const * 100] * len(self.time)
        pl.plot(self.time, Pcnom, '--', color='orange')
        pl.plot(self.time, OtFnom, '--g')
        pl.plot(self.time, OtFx100, 'g', label='O/F Ratio x 100')
        pl.xlabel('Time (s)')
        pl.ylabel('Pressure (psia)')
        pl.grid()
        pl.legend()

        pl.subplot(2, 3, 2)
        pl.title('Engine Thrust')
        pl.plot(self.time, self.thrust, 'hotpink', label='Thrust (lbf)')
        pl.plot(self.time, self.Isp, 'c', label='Specific Impulse (s)')
        thrustnom = [3500] * len(self.time)
        pl.plot(self.time, thrustnom, '--', color='hotpink')
        pl.xlabel('Time (s)')
        pl.ylabel('Thrust (lbf) / Isp (s)')
        pl.grid()
        pl.legend()

        pl.subplot(2, 3, 3)
        pl.title('Added Hydrostatic Pressure')
        pl.plot(self.time, self.P_LOX_hydro, 'b',
                label='LOX Hydrostatic Pressure')
        pl.plot(self.time, self.P_RP1_hydro, 'r',
                label='RP1 Hydrostatic Pressure')
        pl.xlabel('Time (s)')
        pl.ylabel('Pressure (psia)')
        pl.grid()
        pl.legend()

        pl.subplot(2, 3, 4)
        pl.title('Rocket Acceleration')
        Gs = [x / self.g0 for x in self.a]
        pl.plot(self.time, self.a, 'r', label='Rocket Acceleration (m/s^2)')
        pl.plot(self.time, Gs, 'darkred', label='Gs')
        pl.xlabel('Time (s)')
        pl.ylabel('Acceleration (m/s^2)')
        pl.grid()
        pl.legend()

        pl.subplot(2, 3, 5)
        pl.title('Rocket Velocity')
        apo = [0] * len(self.time)
        pl.plot(self.time, self.vel, 'g', label='Rocket Velocity (m/s)')
        pl.plot(self.time, apo, label='Apogee')
        pl.xlabel('Time (s)')
        pl.ylabel('Velocity (m/s)')
        pl.grid()
        pl.legend()

        pl.subplot(2, 3, 6)
        pl.title('Rocket Altitude')
        zkm = [x / 1000 for x in self.z]
        pl.plot(self.time, zkm, 'b', label='Rocket Altitude (km)')
        pl.xlabel('Time (s)')
        pl.ylabel('Altitude (km)')
        pl.grid()
        pl.legend()
        pl.show()

        """

    def plot_animate(self):  # WIP
        return 0


# ----------------------------------------------------------------------------------------------------------------------

    def run(self):
        # Define prop volumes from load profile
        self.set_init_prop()

        # Find CdAs of main prop runlines based on desired pressure inputs, this will become constant as systems mature
        self.set_propCdA()

        # Transient Prepress, good for timing and orifice sizing but otherwise overkill
        self.prepress()

        # Static Fire Burn
        # self.staticfire()

        # Launch Burn and Coasting (Assumes Ballistic)
        self.launch()

        # Ploting or Animated Plotting
        # self.plot()

# ----------------------------------------------------------------------------------------------------------------------


def main():
    test = Halcyon()
    test.run()


if __name__ == '__main__':
    main()
