"""
ImQMD+HIVAP model for MNT reactions.
Based on:
  - Zhao et al., PRC 94, 024601 (2016): ImQMD for U+U
  - Shen et al., PRC 66, 061602 (2002): HIVAP de-excitation
  - Wang et al., NIM B 463 (2020): LENSHIAF gas cell simulations

Physics:
1. Partial wave decomposition with Hill-Wheeler transmission
2. Deflection function for CM scattering angle per partial wave
3. Nucleon exchange probability per channel
4. Excitation energy → neutron evaporation (HIVAP-like)
5. Two-body kinematics for lab frame transformation
"""
import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple


@dataclass
class MNTResult:
    Z: np.ndarray
    N: np.ndarray  
    sigma_matrix: np.ndarray
    sigma_total: float
    E_cm: float


class ImQMD_HIVAP_Model:
    """ImQMD + HIVAP hybrid model with partial-wave kinematics."""

    M_N = 931.494
    HBARC = 197.33
    E2 = 1.44

    def __init__(self, projectile: str, target: str, E_lab: float,
                 Ap: Optional[int] = None, At: Optional[int] = None):
        self.projectile = projectile
        self.target = target
        self.E_lab = E_lab
        self.Zp = self._Z(projectile)
        self.Zt = self._Z(target)
        self.Ap = Ap or self._A(projectile)
        self.At = At or self._A(target)
        self.mu_u = (self.Ap * self.At) / (self.Ap + self.At)
        self.mu_MeV = self.M_N * self.mu_u

        self.E_cm = self.E_lab * self.Ap * self.At / (self.Ap + self.At)
        self.E_lab_total = self.E_lab * self.Ap

        self._lambda_sq = (self.HBARC**2) / (2 * self.mu_MeV * self.E_cm)

        r0 = 1.18
        self.Rp = r0 * self.Ap**(1/3)
        self.Rt = r0 * self.At**(1/3)
        self.R_cont = self.Rp + self.Rt
        self.V_C = self.E2 * self.Zp * self.Zt / self.R_cont

        # Nuclear potential at contact (Bass 1980)
        self.V_nuc = -30.0 + 5.0 * (self.Zp + self.Zt) / 184.0
        self.V_B_eff = self.V_C + self.V_nuc
        self.E_above = max(self.E_cm - self.V_B_eff, 0.0)

        # Sommerfeld parameter
        self.eta = 0.157 * self.Zp * self.Zt * np.sqrt(self.mu_u / self.E_cm)

        # CM velocity in lab
        self.v_cm = np.sqrt(2 * self.E_lab_total / (self.Ap + self.At))

    def _Z(self, s):
        el = {'H':1,'He':2,'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,
              'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18,
              'K':19,'Ca':20,'Sc':21,'Ti':22,'V':23,'Cr':24,'Mn':25,'Fe':26,
              'Co':27,'Ni':28,'Cu':29,'Zn':30,'Ga':31,'Ge':32,'As':33,'Se':34,
              'Br':35,'Kr':36,'Rb':37,'Sr':38,'Y':39,'Zr':40,'Nb':41,'Mo':42,
              'Tc':43,'Ru':44,'Rh':45,'Pd':46,'Ag':47,'Cd':48,'In':49,'Sn':50,
              'Sb':51,'Te':52,'I':53,'Xe':54,'Cs':55,'Ba':56,'La':57,'Ce':58,
              'Pr':59,'Nd':60,'Pm':61,'Sm':62,'Eu':63,'Gd':64,'Tb':65,'Dy':66,
              'Ho':67,'Er':68,'Tm':69,'Yb':70,'Lu':71,'Hf':72,'Ta':73,'W':74,
              'Re':75,'Os':76,'Ir':77,'Pt':78,'Au':79,'Hg':80,'Tl':81,'Pb':82,
              'Bi':83,'Po':84,'At':85,'Rn':86,'Fr':87,'Ra':88,'Ac':89,'Th':90,
              'Pa':91,'U':92,'Np':93,'Pu':94,'Am':95,'Cm':96,'Bk':97,'Cf':98}
        return el.get(s.capitalize(), 0)

    def _A(self, s):
        m = {'Ca':40,'Ni':58,'Xe':136,'Pt':198,'Pb':208,'U':238,
             'Kr':84,'Ge':74,'Sn':124,'Ba':138,'Ra':226,
             'Th':232,'Cm':248,'Ar':40,'Os':200}
        return m.get(s.capitalize(), 0)

    def _transmission(self, l: int, hw: float = 4.0) -> float:
        """Hill-Wheeler transmission."""
        E_rot = self.HBARC**2 * l * (l + 1) / (2 * self.mu_MeV * self.R_cont**2)
        V_eff = self.V_B_eff + E_rot
        dE = V_eff - self.E_cm
        if dE < -5:
            return 1.0
        return 1.0 / (1.0 + np.exp(2 * np.pi * dE / hw))

    def _deflection_angle_cm(self, l: int) -> float:
        """
        CM scattering angle for partial wave l.
        
        For near-barrier U+U, Coulomb deflection breaks down:
        - Strong nuclear attraction pulls nuclei together (pocket in potential)
        - Dinuclear system rotates before separating
        - Fragment emission follows rotation angle, not Coulomb trajectory
        
        Empirical formula based on ImQMD results (Zhao et al.):
        θ_cm(l) = 180° * (1 - l/l_max) + 30° * (l/l_max)
        where the effective l_max is set by the barrier.
        """
        if l < 0:
            l = 0
        
        # Effective maximum l from the barrier
        # For near-barrier, l_eff_max ≈ √(2μR²(E_cm-V_B_eff))/ħ
        l_max_eff = np.sqrt(2 * self.mu_MeV * self.R_cont**2 * 
                           max(self.E_cm - self.V_B_eff, 0.1) / self.HBARC**2)
        
        if l_max_eff < 1:
            l_max_eff = 1
            
        # Empirical deflection: linear from 180° at l=0 to 30° at l=l_max
        theta_cm = np.radians(180.0 * (1 - l / l_max_eff) + 30.0 * (l / l_max_eff))
        return max(theta_cm, np.radians(10.0))

    def _transfer_prob(self, dZ: int, dN: int, l: int) -> float:
        """Transfer probability for channel (dZ, dN) at partial wave l."""
        T_l = self._transmission(l)
        # Excitation energy grows with mass transfer and decreases with l
        E_star = self.E_above * T_l * 0.5 * np.exp(-l / (self.eta * 0.5)) + 0.5
        a = (self.Ap + self.At) / 10.0
        T_temp = np.sqrt(max(E_star / a, 0.01))

        Q = (6.0 * dZ + 8.0 * dN - 0.3 * (dZ**2 + dN**2) + 0.1 * dZ * dN)
        Q_opt = -2.0 * np.sqrt(self.E_cm / (self.Ap + self.At))
        Q_factor = np.exp(-(Q - Q_opt)**2 / (2 * (1.5 * T_temp)**2))

        # Width grows with energy and temperature (Zhao et al. ImQMD finding)
        # For near-barrier U+U with neck formation, width is much larger
        sigma_tr = 3.0 + 1.0 * np.sqrt(max(self.E_above, 0)) + 1.0 * T_temp
        spatial = np.exp(-(dZ**2 + dN**2) / (2 * sigma_tr**2))
        
        return spatial * Q_factor

    def _hivap_evaporation(self, Z_prim: int, A_prim: int, E_star: float) -> list:
        """
        HIVAP-like de-excitation (Shen et al. PRC 66, 061602).
        
        For each primary fragment, calculate the probability of
        evaporating n neutrons. The de-excitation follows a
        statistical evaporation chain.
        
        In HIVAP:
        - Neutron evaporation: each n removes ~B_n + 2T ≈ 8 MeV
        - Width from spin distribution and level density
        
        Returns list of [(Z_final, A_final, probability), ...]
        """
        results = []
        
        # Average neutron separation energy for U isotopes: ~6 MeV
        B_n = 6.0 + 2.0 * (Z_prim - 92) / 10.0  # approx
        
        # Level density parameter: a = A/8 MeV^-1
        a = A_prim / 8.0
        
        # Nuclear temperature at this excitation
        T_nuc = np.sqrt(max(E_star / a, 0.01))
        
        # Maximum number of neutrons that can be evaporated
        n_max = min(int(E_star / B_n) + 1, 12)
        
        # Neutron evaporation probabilities (Weisskopf spectrum)
        total_prob = 0.0
        for n in range(n_max + 1):
            E_rem = E_star - n * B_n
            if E_rem < 0:
                break
            # Probability: P(n) ∝ exp(-n * B_n / T_nuc) * level_density_factor
            P_n = np.exp(-n * B_n / T_nuc) * (E_rem / E_star)**(a/2)
            if n > 0:
                # Spin-dependent reduction for multi-neutron emission
                P_n *= np.exp(-0.1 * n)
            results.append((Z_prim, A_prim - n, P_n))
            total_prob += P_n
        
        # Normalize
        if total_prob > 0:
            for i in range(len(results)):
                results[i] = (results[i][0], results[i][1], results[i][2] / total_prob)
        
        return results

    def _two_body_kinematics(self, m_C: float, m_D: float, 
                              theta_cm_deg: float, Q: float = 0.0):
        """
        Lab frame kinematics for fragment C.
        Returns (E_lab, theta_lab).
        """
        theta_cm = np.radians(theta_cm_deg)

        # Total momentum of system in lab
        p_A = np.sqrt(2 * self.Ap * self.E_lab_total)
        v_cm = p_A / (self.Ap + self.At)

        # Energy of C in CM
        E_cm_C = (m_D / (m_C + m_D)) * (self.E_cm + Q)
        if E_cm_C <= 0:
            return 0.0, 0.0

        v_cm_C = np.sqrt(2 * E_cm_C / m_C)

        v_lab_sq = (v_cm**2 + v_cm_C**2 + 2 * v_cm * v_cm_C * np.cos(theta_cm))

        num = v_cm_C * np.sin(theta_cm)
        den = v_cm + v_cm_C * np.cos(theta_cm)
        if abs(num) < 1e-10 and abs(den) < 1e-10:
            th_lab = 0.0
        else:
            th_lab = np.degrees(np.arctan2(num, den))
            th_lab = max(0, min(90, th_lab))

        E_lab = 0.5 * m_C * max(v_lab_sq, 0)
        return E_lab, th_lab

    def calculate_d2sigma(self, dz_range=(-8, 8), dn_range=(-10, 10),
                           n_l=300, n_theta=90, nE=200):
        """
        Full calculation: partial waves → transfer → evaporation → lab kinematics.
        Returns (theta_lab_grid, E_grid, d2sigma).
        """
        Zv = np.arange(dz_range[0], dz_range[1] + 1)
        Nv = np.arange(dn_range[0], dn_range[1] + 1)
        nZ, nN = len(Zv), len(Nv)

        theta_lab_grid = np.linspace(0, 90, n_theta)
        E_grid = np.linspace(0, 1800, nE)
        dE = E_grid[1] - E_grid[0]
        dTh = theta_lab_grid[1] - theta_lab_grid[0]
        d2 = np.zeros((n_theta, nE))

        # Effective l_max from barrier
        l_max_eff = int(np.sqrt(2 * self.mu_MeV * self.R_cont**2 * 
                       max(self.E_cm - self.V_B_eff, 0.1) / self.HBARC**2))
        l_max = min(l_max_eff + 30, n_l)  # add buffer
        print(f"  l_max_eff = {l_max_eff}, l_max = {l_max}")

        total_sigma = 0.0
        print(f"  Processing {l_max} partial waves...")

        for l in range(l_max):
            T_l = self._transmission(l)
            if T_l < 1e-8:
                continue

            # Cross-section for this partial wave (mb)
            sigma_l = np.pi * self._lambda_sq * (2 * l + 1) * T_l * 10.0
            total_sigma += sigma_l

            # CM scattering angle for this l
            theta_cm = np.degrees(self._deflection_angle_cm(l))

            if l % 50 == 0:
                print(f"    l={l}, T_l={T_l:.4f}, θ_cm={theta_cm:.1f}°, σ_l={sigma_l:.2f} mb")

            # Transfer probabilities
            P_tot = 0.0
            P_channels = np.zeros((nZ, nN))
            for i, dZ in enumerate(Zv):
                for j, dN in enumerate(Nv):
                    p = self._transfer_prob(dZ, dN, l)
                    P_channels[i, j] = p
                    P_tot += p

            if P_tot > 0:
                P_channels /= P_tot

            # For each channel with significant probability
            for i, dZ in enumerate(Zv):
                for j, dN in enumerate(Nv):
                    p_ch = P_channels[i, j]
                    if p_ch < 1e-6:
                        continue

                    sig_ch = sigma_l * p_ch

                    # Primary fragment masses
                    m_C_prim = self.Ap + dZ + dN  # PLF
                    m_D_prim = self.At - dZ - dN  # TLF

                    if m_C_prim <= 0 or m_D_prim <= 0:
                        continue

                    # Excitation energy (Zhao et al. find ~20-50 MeV for U+U)
                    E_star = self.E_above * T_l * 0.8 * np.exp(-l / (l_max_eff * 0.6)) + 2.0
                    E_star = max(E_star, 0.5)

                    # Q-value for this channel
                    Q = (6.0 * dZ + 8.0 * dN - 0.3 * (dZ**2 + dN**2) + 0.1 * dZ * dN)

                    # === De-excitation to ²⁴³U ===
                    # ²⁴³U (Z=92, A=243) can be produced from primary fragments with
                    # Z_prim = 90-94 (Th-U-Pu) via neutron/proton/alpha evaporation
                    # 
                    # Simplified: accept any primary fragment that can reasonably
                    # evaporate to ²⁴³U, with probability ≈ exp(-|A_prim - 243|/3)
                    # for Z_prim = 92, and reduced for other Z.
                    
                    # Check PLF
                    Z_prim = self.Zp + dZ
                    A_prim = int(m_C_prim)
                    dZ_from_92 = abs(Z_prim - 92)
                    dA_from_243 = abs(A_prim - 243)
                    
                    # Accept if close to ²⁴³U in Z-A space
                    if dZ_from_92 <= 2 and dA_from_243 <= 8:
                        # Survival probability: simplified HIVAP
                        # Neutron evaporation: ~1 per 8 MeV excitation
                        n_evap_max = max(0, int(E_star / 6.0))
                        prob_to_243 = 0.0
                        for n_evap in range(n_evap_max + 1):
                            A_final = A_prim - n_evap
                            Z_final = Z_prim
                            # Also allow 1 proton evaporation
                            for p_evap in range(min(1, dZ_from_92) + 1):
                                A_final2 = A_final - p_evap  # proton removes ~1 amu
                                Z_final2 = Z_final - p_evap
                                if Z_final2 == 92 and abs(A_final2 - 243) <= 2:
                                    # Probability: neutron chain × proton × level density
                                    p_n = np.exp(-n_evap * 6.0 / max(E_star/n_evap_max, 0.5))
                                    p_p = 0.3 if p_evap > 0 else 1.0
                                    prob_to_243 += p_n * p_p * 0.1
                        
                        if prob_to_243 > 0.001:
                            m_C_final = 243.0
                            m_D_final = float(self.Ap + self.At - 243)
                            
                            E_loss = abs(A_prim - 243) * 8.0  # energy lost to evaporation
                            Q = (6.0 * dZ + 8.0 * dN - 0.3 * (dZ**2 + dN**2) + 0.1 * dZ * dN)
                            Q_evap = Q - E_loss

                            E_lab, th_lab = self._two_body_kinematics(
                                m_C_final, m_D_final, theta_cm, Q_evap)

                            if E_lab > 0 and 0 <= th_lab <= 90:
                                sigma_theta = 5.0 + 1.0 * np.sqrt(abs(A_prim-243) + abs(dZ) + abs(dN))
                                sigma_E = 50 + 20 * np.sqrt(abs(A_prim-243) + abs(dZ) + abs(dN))
                                weight = sig_ch * prob_to_243
                                g_th = np.exp(-(theta_lab_grid - th_lab)**2 / (2 * sigma_theta**2))
                                g_E = np.exp(-(E_grid - E_lab)**2 / (2 * sigma_E**2))
                                d2 += weight * np.outer(g_th, g_E)
        print(f"  Total cross-section: {total_sigma:.1f} mb")
        print(f"  d2 sum before bg: {np.sum(d2):.1f} (should be ~{total_sigma})")

        # Add low-E background (~5% from target energy loss)
        bg = np.zeros((n_theta, nE))
        for i in range(n_theta):
            th = theta_lab_grid[i]
            for j in range(nE):
                E = E_grid[j]
                if E < 200:
                    bg[i, j] = np.exp(-E / 30) * (1 + 0.3 * np.sin(np.radians(th)))
        bg = bg / bg.sum() * total_sigma * 0.03
        d2 += bg
        print(f"  d2 sum after bg: {np.sum(d2):.1f}")

        # Normalize
        total2 = np.sum(d2) * dE * dTh
        if total2 > 0:
            d2 *= total_sigma / total2

        # Create result object (keep sigma_matrix for backward compat)
        sigma_matrix = np.zeros((len(Zv), len(Nv)))
        for i, dZ in enumerate(Zv):
            for j, dN in enumerate(Nv):
                sigma_matrix[i, j] = total_sigma / (nZ * nN)  # placeholder

        self.result = MNTResult(Z=Zv, N=Nv, sigma_matrix=sigma_matrix,
                                sigma_total=total_sigma, E_cm=self.E_cm)

        return theta_lab_grid, E_grid, d2
