/**
BSD-3-Clause
Copyright 2019 Alliance for Sustainable Energy, LLC
Redistribution and use in source and binary forms, with or without modification, are permitted provided 
that the following conditions are met :
1.	Redistributions of source code must retain the above copyright notice, this list of conditions 
and the following disclaimer.
2.	Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
and the following disclaimer in the documentation and/or other materials provided with the distribution.
3.	Neither the name of the copyright holder nor the names of its contributors may be used to endorse 
or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER, CONTRIBUTORS, UNITED STATES GOVERNMENT OR UNITED STATES 
DEPARTMENT OF ENERGY, NOR ANY OF THEIR EMPLOYEES, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
OR CONSEQUENTIAL DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "csp_solver_gasparticle_tes.h"
#include "csp_solver_util.h"


C_csp_gasparticle_tes::C_csp_gasparticle_tes()
{
	m_vol_tank = m_V_tank_active = m_q_pb_design = m_V_tank_hot_ini = std::numeric_limits<double>::quiet_NaN();

	m_m_dot_tes_dc_max = m_m_dot_tes_ch_max = std::numeric_limits<double>::quiet_NaN();
}

void C_csp_gasparticle_tes::init()
{
    if (!(ms_params.m_ts_hours > 0.0))
        throw(C_csp_exception("The Gen3 gas system must have greater than 0 hours of storage.", "Gas-particle TES initialization"));


	// Declare instance of fluid class for FIELD fluid
	// Set fluid number and copy over fluid matrix if it makes sense
	if( ms_params.m_field_fl != HTFProperties::User_defined && ms_params.m_field_fl < HTFProperties::End_Library_Fluids )
	{
		if( !mc_gas_htfProps.SetFluid(ms_params.m_field_fl) )
		{
			throw(C_csp_exception("Field HTF code is not recognized", "Gas-particle TES Initialization"));
		}
	}
	else if( ms_params.m_field_fl == HTFProperties::User_defined )
	{
		int n_rows = (int)ms_params.m_field_fl_props.nrows();
		int n_cols = (int)ms_params.m_field_fl_props.ncols();
		if( n_rows > 2 && n_cols == 7 )
		{
			if( !mc_gas_htfProps.SetUserDefinedFluid(ms_params.m_field_fl_props) )
			{
				error_msg = util::format(mc_gas_htfProps.UserFluidErrMessage(), n_rows, n_cols);
				throw(C_csp_exception(error_msg, "Two Tank TES Initialization"));
			}
		}
		else
		{
			error_msg = util::format("The user defined field HTF table must contain at least 3 rows and exactly 7 columns. The current table contains %d row(s) and %d column(s)", n_rows, n_cols);
			throw(C_csp_exception(error_msg, "Gas-particle TES Initialization"));
		}
	}
	else
	{
		throw(C_csp_exception("Field HTF code is not recognized", "Gas-particle TES Initialization"));
	}


	// Declare instance of fluid class for STORAGE fluid.
	// Set fluid number and copy over fluid matrix if it makes sense.
	if( ms_params.m_tes_fl != HTFProperties::User_defined && ms_params.m_tes_fl < HTFProperties::End_Library_Fluids )
	{
		if( !mc_particle_htfProps.SetFluid(ms_params.m_tes_fl) )
		{
			throw(C_csp_exception("Storage HTF code is not recognized", "Gas-particle TES Initialization"));
		}	
	}
	else if( ms_params.m_tes_fl == HTFProperties::User_defined )
	{
		int n_rows = (int)ms_params.m_tes_fl_props.nrows();
		int n_cols = (int)ms_params.m_tes_fl_props.ncols();
		if( n_rows > 2 && n_cols == 7 )
		{
			if( !mc_particle_htfProps.SetUserDefinedFluid(ms_params.m_tes_fl_props) )
			{
				error_msg = util::format(mc_particle_htfProps.UserFluidErrMessage(), n_rows, n_cols);
				throw(C_csp_exception(error_msg, "Two Tank TES Initialization"));
			}
		}
		else
		{
			error_msg = util::format("The user defined storage HTF table must contain at least 3 rows and exactly 7 columns. The current table contains %d row(s) and %d column(s)", n_rows, n_cols);
			throw(C_csp_exception(error_msg, "Two Tank TES Initialization"));
		}
	}
	else
	{
		throw(C_csp_exception("Storage HTF code is not recognized", "Two Tank TES Initialization"));
	}

	// Calculate thermal power to PC at design
	m_q_pb_design = ms_params.m_W_dot_pc_design/ms_params.m_eta_pc*1.E6;	//[Wt]

	// Convert parameter units
	ms_params.m_T_cycle_out_des += 273.15;		//[K] convert from C
	ms_params.m_T_receiver_out_des += 273.15;		//[K] convert from C
	ms_params.m_T_tank_hot_ini += 273.15;		//[K] convert from C
	ms_params.m_T_tank_cold_ini += 273.15;		//[K] convert from C


	double Q_tes_des = m_q_pb_design / 1.E6 * ms_params.m_ts_hours;		//[MWt-hr] TES thermal capacity at design

	double q_dot_loss_temp = std::numeric_limits<double>::quiet_NaN();

	gasparticle_tes_sizing(mc_particle_htfProps, Q_tes_des, ms_params.m_T_receiver_out_des - ms_params.m_dt_charge, 
        ms_params.m_T_cycle_out_des + ms_params.m_dt_discharge, ms_params.m_aspect_tank, ms_params.m_tank_pairs, ms_params.m_u_tank,
		m_V_tank_active, m_vol_tank, m_d_silo, m_h_silo, q_dot_loss_temp);

    /*
    initialize the charging heat exchangers. Each receiver is sized to match the power input of the cycle at design, 
    so the heat exchanger duty of each is equal to the design point cycle power.

    Fluid needs to enter the receiver at a temperature of dt_charge (the charging approach temperature) ABOVE the cold
    particle temperature, and the cold particle temperature is equal to T_cycle_out_des + dt_discharge. Therefore, 
    the receiver inlet is at a temperature of T_cycle_out_des + dt_discharge + dt_charge. The charge heat exchanger hot 
    side (gas) needs to restore the flow temperature to this receiver inlet temperature.
    */
    for (size_t i = 0; i < mc_charge_hxs.size(); i++)
        mc_charge_hxs.at(i).init(
                mc_gas_htfProps, mc_particle_htfProps, m_q_pb_design, ms_params.m_dt_charge, 
                ms_params.m_T_receiver_out_des, ms_params.m_T_cycle_out_des + ms_params.m_dt_discharge + ms_params.m_dt_charge
            );
    
    /*
    initialize the discharge heat exchanger. The hot side is the particle flow. 
    */
    mc_discharge_hx.init(
            mc_particle_htfProps, mc_gas_htfProps, m_q_pb_design, ms_params.m_dt_discharge,
            ms_params.m_T_receiver_out_des - ms_params.m_dt_charge, ms_params.m_T_cycle_out_des + ms_params.m_dt_discharge
        );

	// Calculate initial storage values
	double V_inactive = m_vol_tank - m_V_tank_active;
	double V_hot_ini = ms_params.m_f_V_hot_ini*0.01*m_V_tank_active + V_inactive;			//[m^3]
	double V_cold_ini = (1.0 - ms_params.m_f_V_hot_ini*0.01)*m_V_tank_active + V_inactive;	//[m^3]

	double T_hot_ini = ms_params.m_T_tank_hot_ini;		//[K]
	double T_cold_ini = ms_params.m_T_tank_cold_ini;	//[K]

	// Initialize cold and hot tanks
	// Hot tank
	mc_hot_tank.init(mc_particle_htfProps, m_vol_tank, m_h_silo, 0, ms_params.m_u_tank, ms_params.m_tank_pairs, 0, 0, V_hot_ini, T_hot_ini);
	
    // Cold tank
	mc_cold_tank.init(mc_particle_htfProps, m_vol_tank, m_h_silo, 0, ms_params.m_u_tank, ms_params.m_tank_pairs, 0, 0, V_cold_ini, T_cold_ini);

}

bool C_csp_gasparticle_tes::does_tes_exist()
{
	return true;
}

double C_csp_gasparticle_tes::get_hot_temp()
{
	return mc_hot_tank.get_m_T_prev();	//[K]
}

double C_csp_gasparticle_tes::get_cold_temp()
{
	return mc_cold_tank.get_m_T_prev();	//[K]
}

double C_csp_gasparticle_tes::get_initial_charge_energy() 
{
    //MWh
	return m_q_pb_design * ms_params.m_ts_hours * m_V_tank_hot_ini / m_vol_tank *1.e-6;
}

double C_csp_gasparticle_tes::get_min_charge_energy() 
{
    //MWh
    return 0.; //ms_params.m_q_pb_design * ms_params.m_ts_hours * ms_params.m_h_tank_min / ms_params.m_h_tank*1.e-6;
}

double C_csp_gasparticle_tes::get_max_charge_energy() 
{
    //MWh
	//double cp = mc_particle_htfProps.Cp(ms_params.m_T_receiver_out_des);		//[kJ/kg-K] spec heat at average temperature during discharge from hot to cold
 //   double rho = mc_particle_htfProps.dens(ms_params.m_T_receiver_out_des, 1.);

 //   double fadj = (1. - ms_params.m_h_tank_min / ms_params.m_h_tank);

 //   double vol_avail = m_vol_tank * ms_params.m_tank_pairs * fadj;

 //   double e_max = vol_avail * rho * cp * (ms_params.m_T_receiver_out_des - ms_params.m_T_cycle_out_des) / 3.6e6;   //MW-hr

 //   return e_max;
    return m_q_pb_design * ms_params.m_ts_hours / 1.e6;
}

double C_csp_gasparticle_tes::get_degradation_rate()  
{
    //calculates an approximate "average" tank heat loss rate based on some assumptions. Good for simple optimization performance projections.
    double d_tank = sqrt( m_vol_tank / ( (double)ms_params.m_tank_pairs * m_h_silo * 3.14159) );
    double e_loss = ms_params.m_u_tank * 3.14159 * ms_params.m_tank_pairs * d_tank * ( ms_params.m_T_cycle_out_des + ms_params.m_T_receiver_out_des - 576.3 )*1.e-6;  //MJ/s  -- assumes full area for loss, Tamb = 15C
	return e_loss / (m_q_pb_design * ms_params.m_ts_hours * 3600.); //s^-1  -- fraction of heat loss per second based on full charge
}

void C_csp_gasparticle_tes::discharge_avail_est(double T_cold_K, double step_s, double &q_dot_dc_est, double &m_dot_field_est, double &T_hot_field_est) 
{
	double f_storage = 0.0;		// for now, hardcode such that storage always completely discharges

	double m_dot_tank_disch_avail = mc_hot_tank.m_dot_available(f_storage, step_s);	//[kg/s]

    if (m_dot_tank_disch_avail == 0) {
        q_dot_dc_est = 0.;
        m_dot_field_est = 0.;
        T_hot_field_est = std::numeric_limits<double>::quiet_NaN();
        return;
    }

	double T_hot_ini = mc_hot_tank.get_m_T_prev();		//[K]

	double eff, T_cold_tes;
	eff = T_cold_tes = std::numeric_limits<double>::quiet_NaN();
	mc_discharge_hx.hx_discharge_mdot_tes(T_hot_ini, m_dot_tank_disch_avail, T_cold_K, eff, T_cold_tes, T_hot_field_est, q_dot_dc_est, m_dot_field_est);
		
	// If above method fails, it will throw an exception, so if we don't want to break here, need to catch and handle it

	m_m_dot_tes_dc_max = m_dot_field_est;		//[kg/s]
}

void C_csp_gasparticle_tes::charge_avail_est(double T_hot_K, double step_s, double &q_dot_ch_est, double &m_dot_field_est, double &T_cold_field_est)
{
	double f_ch_storage = 0.0;	// for now, hardcode such that storage always completely charges

	double m_dot_tank_charge_avail = mc_cold_tank.m_dot_available(f_ch_storage, step_s);	//[kg/s]

	double T_cold_ini = mc_cold_tank.get_m_T_prev();	//[K]

	double eff, T_hot_tes;
	eff = T_hot_tes = std::numeric_limits<double>::quiet_NaN();

	mc_charge_hxs.front().hx_charge_mdot_tes(T_cold_ini, m_dot_tank_charge_avail, T_hot_K, eff, T_hot_tes, T_cold_field_est, q_dot_ch_est, m_dot_field_est);

	// If above method fails, it will throw an exception, so if we don't want to break here, need to catch and handle it

	m_m_dot_tes_ch_max = m_dot_field_est;		//[kg/s]
}

void C_csp_gasparticle_tes::discharge_full(double timestep /*s*/, double T_amb /*K*/, double T_htf_cold_in /*K*/, double & T_htf_hot_out /*K*/, double & m_dot_htf_out /*kg/s*/, C_csp_tes::S_csp_tes_outputs &outputs)
{
	// This method calculates the timestep-average hot discharge temperature and mass flow rate of the TES system during FULL DISCHARGE.
    // This is out of the field side of the heat exchanger (HX), opposite the tank (or 'TES') side,
    // or if no HX (direct storage), this is equal to the hot tank outlet temperature.

	// Inputs are:
	// 1) Temperature of HTF into TES system. If no heat exchanger, this temperature
	//	   is of the HTF directly entering the cold tank
    
    double q_heater_cold, q_heater_hot, q_dot_loss_cold, q_dot_loss_hot, T_cold_ave, T_hot_ave, m_dot_field, 
        T_field_cold_in, T_field_hot_out, m_dot_tank, T_cold_tank_in;
    q_heater_cold = q_heater_hot = q_dot_loss_cold = q_dot_loss_hot = T_cold_ave = T_hot_ave = m_dot_field 
        = T_field_cold_in = T_field_hot_out = m_dot_tank = T_cold_tank_in = std::numeric_limits<double>::quiet_NaN();

    m_dot_tank = mc_hot_tank.m_dot_available(0.0, timestep);                // [kg/s] maximum tank mass flow for this timestep duration

    mc_hot_tank.energy_balance(timestep, 0.0, m_dot_tank, 0.0, T_amb,       // get average hot tank temperature over timestep
        T_hot_ave, q_heater_hot, q_dot_loss_hot);

    T_field_cold_in = T_htf_cold_in;

    double eff, q_trans;
    eff = q_trans = std::numeric_limits<double>::quiet_NaN();
    mc_discharge_hx.hx_discharge_mdot_tes(T_hot_ave, m_dot_tank, T_field_cold_in,
        eff, T_cold_tank_in, T_field_hot_out, q_trans, m_dot_field);
    
    mc_cold_tank.energy_balance(timestep, m_dot_tank, 0.0, T_cold_tank_in, T_amb,
        T_cold_ave, q_heater_cold, q_dot_loss_cold);

    outputs.m_q_heater = q_heater_hot + q_heater_cold;
    outputs.m_m_dot = m_dot_tank;
    outputs.m_W_dot_rhtf_pump = m_dot_field * ms_params.m_htf_pump_coef / 1.E3 +
        m_dot_tank * ms_params.m_tes_pump_coef / 1.E3;                          //[MWe] Pumping power for Receiver HTF, convert from kW/kg/s*kg/s
    T_htf_hot_out = T_field_hot_out;
    m_dot_htf_out = m_dot_field;
    outputs.m_q_dot_loss = q_dot_loss_hot + q_dot_loss_cold;
    outputs.m_q_dot_ch_from_htf = 0.0;
    outputs.m_T_hot_ave = T_hot_ave;
    outputs.m_T_cold_ave = T_cold_ave;
    outputs.m_T_hot_final = mc_hot_tank.get_m_T_calc();
    outputs.m_T_cold_final = mc_cold_tank.get_m_T_calc();
    
    // Calculate thermal power to HTF
    double T_htf_ave = 0.5*(T_htf_cold_in + T_htf_hot_out);		//[K]
    double cp_htf_ave = mc_gas_htfProps.Cp(T_htf_ave);		//[kJ/kg-K]
    outputs.m_q_dot_dc_to_htf = m_dot_htf_out * cp_htf_ave*(T_htf_hot_out - T_htf_cold_in) / 1000.0;		//[MWt]
}

bool C_csp_gasparticle_tes::discharge(double timestep /*s*/, double T_amb /*K*/, double m_dot_htf_in /*kg/s*/, double T_htf_cold_in /*K*/, double & T_htf_hot_out /*K*/, C_csp_tes::S_csp_tes_outputs &outputs)
{
	// This method calculates the timestep-average hot discharge temperature of the TES system. This is out of the field side of the heat exchanger (HX), opposite the tank (or 'TES') side,
	// or if no HX (direct storage), this is equal to the hot tank outlet temperature.

    // The method returns FALSE if the system output (same as input) mass flow rate is greater than that available

	// Inputs are:
	// 1) Mass flow rate of HTF into TES system (equal to that exiting the system)
	// 2) Temperature of HTF into TES system. If no heat exchanger, this temperature
	//	   is of the HTF directly entering the cold tank

    if (m_dot_htf_in > m_m_dot_tes_dc_max)      // mass flow in = mass flow out
    {
        outputs.m_q_heater = std::numeric_limits<double>::quiet_NaN();
        outputs.m_m_dot = std::numeric_limits<double>::quiet_NaN();
        outputs.m_W_dot_rhtf_pump = std::numeric_limits<double>::quiet_NaN();
        outputs.m_q_dot_loss = std::numeric_limits<double>::quiet_NaN();
        outputs.m_q_dot_dc_to_htf = std::numeric_limits<double>::quiet_NaN();
        outputs.m_q_dot_ch_from_htf = std::numeric_limits<double>::quiet_NaN();
        outputs.m_T_hot_ave = std::numeric_limits<double>::quiet_NaN();
        outputs.m_T_cold_ave = std::numeric_limits<double>::quiet_NaN();
        outputs.m_T_hot_final = std::numeric_limits<double>::quiet_NaN();
        outputs.m_T_cold_final = std::numeric_limits<double>::quiet_NaN();

        return false;
    }

	double q_heater_cold, q_heater_hot, q_dot_loss_cold, q_dot_loss_hot, T_cold_ave, T_hot_ave, m_dot_field, T_field_cold_in, T_field_hot_out, m_dot_tank, T_cold_tank_in;
    q_heater_cold = q_heater_hot = q_dot_loss_cold = q_dot_loss_hot = T_cold_ave = T_hot_ave = m_dot_field = T_field_cold_in = T_field_hot_out = m_dot_tank = T_cold_tank_in = std::numeric_limits<double>::quiet_NaN();

	// Iterate between field htf - hx - and storage
    m_dot_field = m_dot_htf_in;
    T_field_cold_in = T_htf_cold_in;
    C_MEQ_indirect_tes_discharge c_tes_disch(this, timestep, T_amb, T_field_cold_in, m_dot_field);
    C_monotonic_eq_solver c_tes_disch_solver(c_tes_disch);

    // Set up solver for tank mass flow
    double m_dot_tank_lower = 0;
    double m_dot_tank_upper = mc_hot_tank.m_dot_available(0.0, timestep);
    c_tes_disch_solver.settings(1.E-3, 50, m_dot_tank_lower, m_dot_tank_upper, false);

    // Guess tank mass flows
    double T_tank_guess_1, m_dot_tank_guess_1, T_tank_guess_2, m_dot_tank_guess_2;
    T_tank_guess_1 = m_dot_tank_guess_1 = T_tank_guess_2 = m_dot_tank_guess_2 = std::numeric_limits<double>::quiet_NaN();
    double eff_guess, T_hot_field_guess, T_cold_tes_guess, q_trans_guess;
    eff_guess = T_hot_field_guess = T_cold_tes_guess = q_trans_guess = std::numeric_limits<double>::quiet_NaN();

    T_tank_guess_1 = mc_hot_tank.get_m_T_prev();                                // the highest possible tank temp and thus highest resulting mass flow
    mc_discharge_hx.hx_discharge_mdot_field(T_field_cold_in, m_dot_field, T_tank_guess_1,
        eff_guess, T_hot_field_guess, T_cold_tes_guess, q_trans_guess, m_dot_tank_guess_1);

    double T_s, mCp, Q_u, L_s, UA;
    T_s = mc_hot_tank.get_m_T_prev();
    mCp = mc_hot_tank.calc_mass_at_prev() * mc_hot_tank.calc_cp_at_prev();  // [J/K]
    Q_u = 0.;   // [W]
    //L_s = mc_hot_tank.m_dot_available(0.0, timestep) * mc_hot_tank.calc_enth_at_prev();       // maybe use this instead?
    L_s = m_dot_tank_guess_1 * mc_hot_tank.calc_enth_at_prev();  // [W]  using highest mass flow, calculated from the first guess
    UA = mc_hot_tank.get_m_UA();
    T_tank_guess_2 = T_s + timestep / mCp * (Q_u - L_s - UA * (T_s - T_amb));   // tank energy balance, giving the lowest estimated tank temp due to L_s choice
    mc_discharge_hx.hx_discharge_mdot_field(T_field_cold_in, m_dot_field, T_tank_guess_2,
        eff_guess, T_hot_field_guess, T_cold_tes_guess, q_trans_guess, m_dot_tank_guess_2);

    // Constrain guesses to within limits
    m_dot_tank_guess_1 = fmax(m_dot_tank_lower, fmin(m_dot_tank_guess_1, m_dot_tank_upper));
    m_dot_tank_guess_2 = fmax(m_dot_tank_lower, fmin(m_dot_tank_guess_2, m_dot_tank_upper));

    bool m_dot_solved = false;
    if (m_dot_tank_guess_1 == m_dot_tank_guess_2) {
        // First try if guess solves
        double m_dot_bal = std::numeric_limits<double>::quiet_NaN();
        int m_dot_bal_code = c_tes_disch_solver.test_member_function(m_dot_tank_guess_1, &m_dot_bal);

        if (m_dot_bal < 1.E-3) {
            m_dot_tank = m_dot_tank_guess_1;
            m_dot_solved = true;
        }
        else {
            // Adjust guesses so they're different
            if (m_dot_tank_guess_1 == m_dot_tank_upper) {
                m_dot_tank_guess_2 -= 0.01*(m_dot_tank_upper - m_dot_tank_lower);
            }
            else {
                m_dot_tank_guess_1 = fmin(m_dot_tank_upper, m_dot_tank_guess_1 + 0.01*(m_dot_tank_upper - m_dot_tank_lower));
            }
        }
    }
    else if (m_dot_tank_guess_1 == 0 || m_dot_tank_guess_2 == 0) {
        // Try a 0 guess
        double m_dot_bal = std::numeric_limits<double>::quiet_NaN();
        int m_dot_bal_code = c_tes_disch_solver.test_member_function(0., &m_dot_bal);

        if (m_dot_bal < 1.E-3) {
            m_dot_tank = 0.;
            m_dot_solved = true;
        }
        else {
            // Adjust 0 guess to avoid divide by 0 errors
            m_dot_tank_guess_2 = fmax(m_dot_tank_guess_1, m_dot_tank_guess_2);
            m_dot_tank_guess_1 = 1.e-3;
        }
    }

    if (!m_dot_solved) {
        // Solve for required tank mass flow
        double tol_solved;
        tol_solved = std::numeric_limits<double>::quiet_NaN();
        int iter_solved = -1;

        try
        {
            int m_dot_bal_code = c_tes_disch_solver.solve(m_dot_tank_guess_1, m_dot_tank_guess_2, -1.E-3, m_dot_tank, tol_solved, iter_solved);
        }
        catch (C_csp_exception)
        {
            throw(C_csp_exception("Failed to find a solution for the hot tank mass flow"));
        }

        if (std::isnan(m_dot_tank)) {
            throw(C_csp_exception("Failed to converge on a valid tank mass flow."));
        }
    }

    // The other needed outputs are not all saved in member variables so recalculate here
    mc_hot_tank.energy_balance(timestep, 0.0, m_dot_tank, 0.0, T_amb, T_hot_ave, q_heater_hot, q_dot_loss_hot);

    double eff, q_trans;
    eff = q_trans = std::numeric_limits<double>::quiet_NaN();
    mc_discharge_hx.hx_discharge_mdot_field(T_field_cold_in, m_dot_field, T_hot_ave,
        eff, T_field_hot_out, T_cold_tank_in, q_trans, m_dot_tank);

        mc_cold_tank.energy_balance(timestep, m_dot_tank, 0.0, T_cold_tank_in, T_amb, T_cold_ave, q_heater_cold, q_dot_loss_cold);

	outputs.m_q_heater = q_heater_cold + q_heater_hot;			//[MWt]
    outputs.m_m_dot = m_dot_tank;                               //[kg/s]
    outputs.m_W_dot_rhtf_pump = m_dot_field * ms_params.m_htf_pump_coef / 1.E3 +
        m_dot_tank * ms_params.m_tes_pump_coef / 1.E3;      //[MWe] Pumping power for Receiver HTF, convert from kW/kg/s*kg/s
    T_htf_hot_out = T_field_hot_out;
	outputs.m_q_dot_loss = q_dot_loss_cold + q_dot_loss_hot;	//[MWt]
	outputs.m_q_dot_ch_from_htf = 0.0;		                    //[MWt]
	outputs.m_T_hot_ave = T_hot_ave;						    //[K]
	outputs.m_T_cold_ave = T_cold_ave;							//[K]
	outputs.m_T_hot_final = mc_hot_tank.get_m_T_calc();			//[K]
	outputs.m_T_cold_final = mc_cold_tank.get_m_T_calc();		//[K]

	// Calculate thermal power to HTF
	double T_htf_ave = 0.5*(T_htf_cold_in + T_htf_hot_out);		//[K]
	double cp_htf_ave = mc_gas_htfProps.Cp(T_htf_ave);		//[kJ/kg-K]
	outputs.m_q_dot_dc_to_htf = m_dot_htf_in*cp_htf_ave*(T_htf_hot_out - T_htf_cold_in)/1000.0;		//[MWt]

	return true;
}

bool C_csp_gasparticle_tes::charge(C_heat_exchanger& charge_hx, double timestep /*s*/, double T_amb /*K*/, double m_dot_htf_in /*kg/s*/, double T_htf_hot_in /*K*/, double & T_htf_cold_out /*K*/, C_csp_tes::S_csp_tes_outputs &outputs)
{
	// This method calculates the timestep-average cold charge return temperature of the TES system. This is out of the field side of the heat exchanger (HX), opposite the tank (or 'TES') side,
    // or if no HX (direct storage), this is equal to the cold tank outlet temperature.
	
	// The method returns FALSE if the system input mass flow rate is greater than the allowable charge 

	// Inputs are:
	// 1) Mass flow rate of HTF into TES system (equal to that exiting the system)
	// 2) Temperature of HTF into TES system. If no heat exchanger, this temperature
	//	   is of the HTF directly entering the hot tank

    if (m_dot_htf_in > m_m_dot_tes_ch_max)
    {
        outputs.m_q_heater = std::numeric_limits<double>::quiet_NaN();
        outputs.m_m_dot = std::numeric_limits<double>::quiet_NaN();
        outputs.m_W_dot_rhtf_pump = std::numeric_limits<double>::quiet_NaN();
        outputs.m_q_dot_loss = std::numeric_limits<double>::quiet_NaN();
        outputs.m_q_dot_dc_to_htf = std::numeric_limits<double>::quiet_NaN();
        outputs.m_q_dot_ch_from_htf = std::numeric_limits<double>::quiet_NaN();
        outputs.m_T_hot_ave = std::numeric_limits<double>::quiet_NaN();
        outputs.m_T_cold_ave = std::numeric_limits<double>::quiet_NaN();
        outputs.m_T_hot_final = std::numeric_limits<double>::quiet_NaN();
        outputs.m_T_cold_final = std::numeric_limits<double>::quiet_NaN();

        return false;
    }

	double q_heater_cold, q_heater_hot, q_dot_loss_cold, q_dot_loss_hot, T_cold_ave, T_hot_ave, m_dot_field, T_field_hot_in, T_field_cold_out, m_dot_tank, T_hot_tank_in;
	q_heater_cold = q_heater_hot = q_dot_loss_cold = q_dot_loss_hot = T_cold_ave = T_hot_ave = m_dot_field = T_field_hot_in = T_field_cold_out = m_dot_tank = T_hot_tank_in = std::numeric_limits<double>::quiet_NaN();

	// If no heat exchanger, no iteration is required between the heat exchanger and storage tank models
	// Iterate between field htf - hx - and storage	
    m_dot_field = m_dot_htf_in;
    T_field_hot_in = T_htf_hot_in;
    C_MEQ_indirect_tes_charge c_tes_chrg(this, timestep, T_amb, T_field_hot_in, m_dot_field);
    C_monotonic_eq_solver c_tes_chrg_solver(c_tes_chrg);

    // Set up solver for tank mass flow
    double m_dot_tank_lower = 0;
    double m_dot_tank_upper = mc_cold_tank.m_dot_available(0.0, timestep);
    c_tes_chrg_solver.settings(1.E-3, 50, m_dot_tank_lower, m_dot_tank_upper, false);

    // Guess tank mass flows
    double T_tank_guess_1, m_dot_tank_guess_1, T_tank_guess_2, m_dot_tank_guess_2;
    T_tank_guess_1 = m_dot_tank_guess_1 = T_tank_guess_2 = m_dot_tank_guess_2 = std::numeric_limits<double>::quiet_NaN();
    double eff_guess, T_cold_field_guess, T_hot_tes_guess, q_trans_guess;
    eff_guess = T_cold_field_guess = T_hot_tes_guess = q_trans_guess = std::numeric_limits<double>::quiet_NaN();
    
    T_tank_guess_1 = mc_cold_tank.get_m_T_prev();                                // the highest possible tank temp and thus highest resulting mass flow
    charge_hx.hx_charge_mdot_field(T_field_hot_in, m_dot_field, T_tank_guess_1,
        eff_guess, T_cold_field_guess, T_hot_tes_guess, q_trans_guess, m_dot_tank_guess_1);

    double T_s, mCp, Q_u, L_s, UA;
    T_s = mc_cold_tank.get_m_T_prev();
    mCp = mc_cold_tank.calc_mass_at_prev() * mc_cold_tank.calc_cp_at_prev();  // [J/K]
    Q_u = 0.;   // [W]
    //L_s = mc_cold_tank.m_dot_available(0.0, timestep) * mc_cold_tank.calc_enth_at_prev();       // maybe use this instead?
    L_s = m_dot_tank_guess_1 * mc_cold_tank.calc_enth_at_prev();  // [W]  using highest mass flow, calculated from the first guess
    UA = mc_cold_tank.get_m_UA();
    T_tank_guess_2 = T_s + timestep / mCp * (Q_u - L_s - UA * (T_s - T_amb));   // tank energy balance, giving the lowest estimated tank temp due to L_s choice
    charge_hx.hx_charge_mdot_field(T_field_hot_in, m_dot_field, T_tank_guess_2,
        eff_guess, T_cold_field_guess, T_hot_tes_guess, q_trans_guess, m_dot_tank_guess_2);

    // Constrain guesses to within limits
    m_dot_tank_guess_1 = fmax(m_dot_tank_lower, fmin(m_dot_tank_guess_1, m_dot_tank_upper));
    m_dot_tank_guess_2 = fmax(m_dot_tank_lower, fmin(m_dot_tank_guess_2, m_dot_tank_upper));

    bool m_dot_solved = false;
    if (m_dot_tank_guess_1 == m_dot_tank_guess_2) {
        // First try if guess solves
        double m_dot_bal = std::numeric_limits<double>::quiet_NaN();
        int m_dot_bal_code = c_tes_chrg_solver.test_member_function(m_dot_tank_guess_1, &m_dot_bal);

        if (m_dot_bal < 1.E-3) {
            m_dot_tank = m_dot_tank_guess_1;
            m_dot_solved = true;
        }
        else {
            // Adjust guesses so they're different
            if (m_dot_tank_guess_1 == m_dot_tank_upper) {
                m_dot_tank_guess_2 -= 0.01*(m_dot_tank_upper - m_dot_tank_lower);
            }
            else {
                m_dot_tank_guess_1 = fmin(m_dot_tank_upper, m_dot_tank_guess_1 + 0.01*(m_dot_tank_upper - m_dot_tank_lower));
            }
        }
    }
    else if (m_dot_tank_guess_1 == 0 || m_dot_tank_guess_2 == 0) {
        // Try a 0 guess
        double m_dot_bal = std::numeric_limits<double>::quiet_NaN();
        int m_dot_bal_code = c_tes_chrg_solver.test_member_function(0., &m_dot_bal);

        if (m_dot_bal < 1.E-3) {
            m_dot_tank = 0.;
            m_dot_solved = true;
        }
        else {
            // Adjust 0 guess to avoid divide by 0 errors
            m_dot_tank_guess_2 = fmax(m_dot_tank_guess_1, m_dot_tank_guess_2);
            m_dot_tank_guess_1 = 1.e-3;
        }
    }

    if (!m_dot_solved) {
        // Solve for required tank mass flow
        double tol_solved;
        tol_solved = std::numeric_limits<double>::quiet_NaN();
        int iter_solved = -1;

        try
        {
            int m_dot_bal_code = c_tes_chrg_solver.solve(m_dot_tank_guess_1, m_dot_tank_guess_2, -1.E-3, m_dot_tank, tol_solved, iter_solved);
        }
        catch (C_csp_exception)
        {
            throw(C_csp_exception("Failed to find a solution for the cold tank mass flow"));
        }

        if (std::isnan(m_dot_tank)) {
            throw(C_csp_exception("Failed to converge on a valid tank mass flow."));
        }
    }

    // The other needed outputs are not all saved in member variables so recalculate here
    mc_cold_tank.energy_balance(timestep, 0.0, m_dot_tank, 0.0, T_amb, T_cold_ave, q_heater_cold, q_dot_loss_cold);

    double eff, q_trans;
    eff = q_trans = std::numeric_limits<double>::quiet_NaN();
    charge_hx.hx_charge_mdot_field(T_field_hot_in, m_dot_field, T_cold_ave,
        eff, T_field_cold_out, T_hot_tank_in, q_trans, m_dot_tank);

    mc_hot_tank.energy_balance(timestep, m_dot_tank, 0.0, T_hot_tank_in, T_amb, T_hot_ave, q_heater_hot, q_dot_loss_hot);

	outputs.m_q_heater = q_heater_cold + q_heater_hot;			//[MWt]
    outputs.m_m_dot = m_dot_tank;                               //[kg/s]
    outputs.m_W_dot_rhtf_pump = m_dot_field * ms_params.m_htf_pump_coef / 1.E3 +
        m_dot_tank * ms_params.m_tes_pump_coef / 1.E3;      //[MWe] Pumping power for Receiver HTF, convert from kW/kg/s*kg/s
    T_htf_cold_out = T_field_cold_out;
    outputs.m_q_dot_loss = q_dot_loss_cold + q_dot_loss_hot;	//[MWt]
	outputs.m_q_dot_dc_to_htf = 0.0;							//[MWt]
	outputs.m_T_hot_ave = T_hot_ave;							//[K]
	outputs.m_T_cold_ave = T_cold_ave;  						//[K]
	outputs.m_T_hot_final = mc_hot_tank.get_m_T_calc();			//[K]
	outputs.m_T_cold_final = mc_cold_tank.get_m_T_calc();		//[K]

	// Calculate thermal power to HTF
	double T_htf_ave = 0.5*(T_htf_hot_in + T_htf_cold_out);		//[K]
	double cp_htf_ave = mc_gas_htfProps.Cp(T_htf_ave);		//[kJ/kg-K]
	outputs.m_q_dot_ch_from_htf = m_dot_htf_in*cp_htf_ave*(T_htf_hot_in - T_htf_cold_out)/1000.0;		//[MWt]

	return true;
}

void C_csp_gasparticle_tes::charge_full(double timestep /*s*/, double T_amb /*K*/, double T_htf_hot_in /*K*/, double & T_htf_cold_out /*K*/, double & m_dot_htf_out /*kg/s*/, C_csp_tes::S_csp_tes_outputs &outputs)
{
	// This method calculates the timestep-average cold charge return temperature and mass flow rate of the TES system during during FULL CHARGE.
    // This is out of the field side of the heat heat exchanger (HX), opposite the tank (or 'TES') side,
	// or if no HX (direct storage), this is equal to the cold tank outlet temperature.

	// Inputs are:
	// 1) Temperature of HTF into TES system. If no heat exchanger, this temperature
	//	   is of the HTF directly entering the hot tank

    double q_heater_cold, q_heater_hot, q_dot_loss_cold, q_dot_loss_hot, T_cold_ave, T_hot_ave, m_dot_field, T_field_cold_in, T_field_hot_out, m_dot_tank, T_hot_tank_in;
    q_heater_cold = q_heater_hot = q_dot_loss_cold = q_dot_loss_hot = T_cold_ave = T_hot_ave = m_dot_field = T_field_cold_in = T_field_hot_out = m_dot_tank = T_hot_tank_in = std::numeric_limits<double>::quiet_NaN();

    m_dot_tank = mc_cold_tank.m_dot_available(0.0, timestep);               // [kg/s] maximum tank mass flow for this timestep duration

    mc_cold_tank.energy_balance(timestep, 0.0, m_dot_tank, 0.0, T_amb,      // get average hot tank temperature over timestep
        T_cold_ave, q_heater_cold, q_dot_loss_cold);

    // Get hot tank inlet temperature
    T_field_hot_out = T_htf_hot_in;

    double eff, q_trans;
    eff = q_trans = std::numeric_limits<double>::quiet_NaN();
    mc_charge_hxs.front().hx_charge_mdot_tes(T_cold_ave, m_dot_tank / (double)mc_charge_hxs.size(), T_field_hot_out,
        eff, T_hot_tank_in, T_field_cold_in, q_trans, m_dot_field);

    mc_hot_tank.energy_balance(timestep, m_dot_tank, 0.0, T_hot_tank_in, T_amb,
        T_hot_ave, q_heater_hot, q_dot_loss_hot);

    outputs.m_q_heater = q_heater_hot + q_heater_cold;
    outputs.m_m_dot = m_dot_tank;
    outputs.m_W_dot_rhtf_pump = m_dot_field * ms_params.m_htf_pump_coef / 1.E3 +
        m_dot_tank * ms_params.m_tes_pump_coef / 1.E3;                          //[MWe] Pumping power for Receiver HTF, convert from kW/kg/s*kg/s
    T_htf_cold_out = T_field_cold_in;
    m_dot_htf_out = m_dot_field;
    outputs.m_q_dot_loss = q_dot_loss_hot + q_dot_loss_cold;
    outputs.m_q_dot_dc_to_htf = 0.0;
    outputs.m_T_hot_ave = T_hot_ave;
    outputs.m_T_cold_ave = T_cold_ave;
    outputs.m_T_hot_final = mc_hot_tank.get_m_T_calc();
    outputs.m_T_cold_final = mc_cold_tank.get_m_T_calc();

    // Calculate thermal power to HTF
    double T_htf_ave = 0.5*(T_htf_hot_in + T_htf_cold_out);		//[K]
    double cp_htf_ave = mc_gas_htfProps.Cp(T_htf_ave);		//[kJ/kg-K]
    outputs.m_q_dot_ch_from_htf = m_dot_htf_out * cp_htf_ave*(T_htf_hot_in - T_htf_cold_out) / 1000.0;		//[MWt]
}

void C_csp_gasparticle_tes::idle(double timestep, double T_amb, C_csp_tes::S_csp_tes_outputs &outputs)
{
	double T_hot_ave, q_hot_heater, q_dot_hot_loss;
	T_hot_ave = q_hot_heater = q_dot_hot_loss = std::numeric_limits<double>::quiet_NaN();

	mc_hot_tank.energy_balance(timestep, 0.0, 0.0, 0.0, T_amb, T_hot_ave, q_hot_heater, q_dot_hot_loss);

	double T_cold_ave, q_cold_heater, q_dot_cold_loss;
	T_cold_ave = q_cold_heater = q_dot_cold_loss = std::numeric_limits<double>::quiet_NaN();

	mc_cold_tank.energy_balance(timestep, 0.0, 0.0, 0.0, T_amb, T_cold_ave, q_cold_heater, q_dot_cold_loss);

	outputs.m_q_heater = q_cold_heater + q_hot_heater;			//[MJ]
	outputs.m_W_dot_rhtf_pump = 0.0;							//[MWe]
	outputs.m_q_dot_loss = q_dot_cold_loss + q_dot_hot_loss;	//[MW]
	
	outputs.m_T_hot_ave = T_hot_ave;							//[K]
	outputs.m_T_cold_ave = T_cold_ave;							//[K]
	outputs.m_T_hot_final = mc_hot_tank.get_m_T_calc();			//[K]
	outputs.m_T_cold_final = mc_cold_tank.get_m_T_calc();		//[K]

	outputs.m_q_dot_ch_from_htf = 0.0;		//[MWt]
	outputs.m_q_dot_dc_to_htf = 0.0;		//[MWt]
}

void C_csp_gasparticle_tes::converged()
{
	mc_cold_tank.converged();
	mc_hot_tank.converged();
    for(size_t i=0; i<mc_charge_hxs.size(); i++)
        mc_charge_hxs.at(i).converged();
	
	// The max charge and discharge flow rates should be set at the beginning of each timestep
	//   during the q_dot_xx_avail_est calls
	m_m_dot_tes_dc_max = m_m_dot_tes_ch_max = std::numeric_limits<double>::quiet_NaN();
}

void gasparticle_tes_sizing(HTFProperties &tes_htf_props, double Q_tes_des /*MWt-hr*/, double T_tes_hot /*K*/,
	double T_tes_cold /*K*/, double aspect_silo /*m*/, int silo_pairs /*-*/, double u_silo /*W/m^2-K*/,
	double & vol_one_temp_avail /*m3*/, double & vol_one_temp_total /*m3*/, double & d_silo /*m*/, double & h_silo /*m*/,
	double & q_dot_loss_des /*MWt*/)
{
	double T_tes_ave = 0.5*(T_tes_hot + T_tes_cold);		//[K]
	
	double rho_ave = tes_htf_props.dens(T_tes_ave, 1.0);		//[kg/m^3] Density at average temperature
	double cp_ave = tes_htf_props.Cp(T_tes_ave);				//[kJ/kg-K] Specific heat at average temperature

	// Volume required to supply design hours of thermal energy storage
		//[m^3] = [MJ/s-hr] * [sec]/[hr] = [MJ] / (kg/m^3 * MJ/kg-K * K 
	vol_one_temp_avail = Q_tes_des*3600.0 / (rho_ave * cp_ave / 1000.0 * (T_tes_hot - T_tes_cold));

	// Additional volume necessary due to minimum silo limits
	vol_one_temp_total = vol_one_temp_avail / (1.0 - 0.1 / h_silo);	//[m^3]  Put in a guess that the tank is 10% oversized for the particle volume

    // solve for tank dimensions given volume and aspect
    d_silo = pow(4.* vol_one_temp_total / (CSP::pi * aspect_silo * silo_pairs), 1. / 3.);
    h_silo = d_silo * aspect_silo;


	double A_cs = vol_one_temp_total / (h_silo*silo_pairs);		//[m^2] Cross-sectional area of a single silo

	d_silo = pow(A_cs / CSP::pi, 0.5)*2.0;			//[m] Diameter of a single silo

	double UA_silo = u_silo*(A_cs + CSP::pi*d_silo*h_silo)*silo_pairs;		//[W/K]

	q_dot_loss_des = UA_silo*(T_tes_ave - 15.0)*1.E-6;	//[MWt]
		
}

int C_csp_gasparticle_tes::C_MEQ_indirect_tes_discharge::operator()(double m_dot_tank /*kg/s*/, double *m_dot_bal /*-*/)
{
    // Call energy balance on hot tank discharge using a guess for storage-side mass flow to get average tank outlet temperature over timestep
    double T_hot_ave, q_heater_hot, q_dot_loss_hot;
    T_hot_ave = q_heater_hot = q_dot_loss_hot = std::numeric_limits<double>::quiet_NaN();
    mpc_csp_gasparticle_tes->mc_hot_tank.energy_balance(m_timestep, 0.0, m_dot_tank, 0.0, m_T_amb, T_hot_ave, q_heater_hot, q_dot_loss_hot);

    // Use average tank outlet temperature to call energy balance on heat exchanger to get storage-side mass flow
    double eff, T_hot_field, T_cold_tes, q_trans, m_dot_tank_solved;
    eff = T_hot_field = T_cold_tes = q_trans = std::numeric_limits<double>::quiet_NaN();
    mpc_csp_gasparticle_tes->mc_discharge_hx.hx_discharge_mdot_field(m_T_cold_field, m_m_dot_field, T_hot_ave,
        eff, T_hot_field, T_cold_tes, q_trans, m_dot_tank_solved);

    if (m_dot_tank != 0.) {
        *m_dot_bal = (m_dot_tank_solved - m_dot_tank) / m_dot_tank;			//[-]
    }
    else {
        *m_dot_bal = m_dot_tank_solved - m_dot_tank;            			//[kg/s]  return absolute difference if 0
    }

    return 0;
}

int C_csp_gasparticle_tes::C_MEQ_indirect_tes_charge::operator()(double m_dot_tank /*kg/s*/, double *m_dot_bal /*-*/)
{
    // Call energy balance on cold tank discharge using a guess for storage-side mass flow to get average tank outlet temperature over timestep
    double T_cold_ave, q_heater_cold, q_dot_loss_cold;
    T_cold_ave = q_heater_cold = q_dot_loss_cold = std::numeric_limits<double>::quiet_NaN();
    mpc_csp_gasparticle_tes->mc_cold_tank.energy_balance(m_timestep, 0.0, m_dot_tank, 0.0, m_T_amb, T_cold_ave, q_heater_cold, q_dot_loss_cold);

    // Use average tank outlet temperature to call energy balance on heat exchanger to get storage-side mass flow
    double eff, T_cold_field, T_hot_tes, q_trans, m_dot_tank_solved;
    eff = T_cold_field = T_hot_tes = q_trans = std::numeric_limits<double>::quiet_NaN();
    
    mpc_csp_gasparticle_tes->mc_charge_hxs.front().hx_charge_mdot_field(m_T_hot_field, m_m_dot_field, T_cold_ave,
        eff, T_cold_field, T_hot_tes, q_trans, m_dot_tank_solved);

    if (m_dot_tank != 0.) {
        *m_dot_bal = (m_dot_tank_solved - m_dot_tank) / m_dot_tank;			//[-]
    }
    else {
        *m_dot_bal = m_dot_tank_solved - m_dot_tank;            			//[kg/s]  return absolute difference if 0
    }

    return 0;
}

