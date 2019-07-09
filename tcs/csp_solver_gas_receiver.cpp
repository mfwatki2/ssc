/*******************************************************************************************************
*  Copyright 2017 Alliance for Sustainable Energy, LLC
*
*  NOTICE: This software was developed at least in part by Alliance for Sustainable Energy, LLC
*  (“Alliance”) under Contract No. DE-AC36-08GO28308 with the U.S. Department of Energy and the U.S.
*  The Government retains for itself and others acting on its behalf a nonexclusive, paid-up,
*  irrevocable worldwide license in the software to reproduce, prepare derivative works, distribute
*  copies to the public, perform publicly and display publicly, and to permit others to do so.
*
*  Redistribution and use in source and binary forms, with or without modification, are permitted
*  provided that the following conditions are met:
*
*  1. Redistributions of source code must retain the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer.
*
*  2. Redistributions in binary form must reproduce the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer in the documentation and/or
*  other materials provided with the distribution.
*
*  3. The entire corresponding source code of any redistribution, with or without modification, by a
*  research entity, including but not limited to any contracting manager/operator of a United States
*  National Laboratory, any institution of higher learning, and any non-profit organization, must be
*  made publicly available under this license for as long as the redistribution is made available by
*  the research entity.
*
*  4. Redistribution of this software, without modification, must refer to the software by the same
*  designation. Redistribution of a modified version of this software (i) may not refer to the modified
*  version by the same designation, or by any confusingly similar designation, and (ii) must refer to
*  the underlying software originally provided by Alliance as “System Advisor Model” or “SAM”. Except
*  to comply with the foregoing, the terms “System Advisor Model”, “SAM”, or any confusingly similar
*  designation may not be used to refer to any modified version of this software or any modified
*  version of the underlying software originally provided by Alliance without the prior written consent
*  of Alliance.
*
*  5. The name of the copyright holder, contributors, the United States Government, the United States
*  Department of Energy, or any of their employees may not be used to endorse or promote products
*  derived from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
*  IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
*  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER,
*  CONTRIBUTORS, UNITED STATES GOVERNMENT OR UNITED STATES DEPARTMENT OF ENERGY, NOR ANY OF THEIR
*  EMPLOYEES, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
*  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
*  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
*  IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
*  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************************************/

#include "csp_solver_gas_receiver.h"
#include "csp_solver_gas.h"
#include "sam_csp_util.h"


/*
	NodeNodeNode
	[..][..][..]
<<	------------      --^--
	||||||||||||        |
	||||||||||||        |
	||||||||||||        |
>>	============	  h_rec
	||||||||||||        |
	||||||||||||        |
	||||||||||||        |
<<	------------       _v_

Receiver consists of a bunch of tubes branching off of a centerline header. Tubes go either vertically up or down and 
flow is collected at top and bottom of the receiver.

* A node consists of a set of parallel tubes
* Nodes are symmetric about centerline
	* Only model one half of the panel and assume symmetry
	* Mass flow in is divided by half for top/bottom split
* Nodes have shared inlet/outlet pressure
* The number of nodes is equal to m_n_flux_x

*/

C_gen3gas_receiver::C_gen3gas_receiver()
{
	m_w_rec = std::numeric_limits<double>::quiet_NaN();
	m_h_rec = std::numeric_limits<double>::quiet_NaN();
	m_h_tower = std::numeric_limits<double>::quiet_NaN();
	m_od_tube = std::numeric_limits<double>::quiet_NaN();
	m_th_tube = std::numeric_limits<double>::quiet_NaN();
	m_epsilon = std::numeric_limits<double>::quiet_NaN();
	m_hl_ffact = std::numeric_limits<double>::quiet_NaN();
	m_T_hot_des = std::numeric_limits<double>::quiet_NaN();
	m_T_cold_des = std::numeric_limits<double>::quiet_NaN();
	m_f_rec_min = std::numeric_limits<double>::quiet_NaN();
	m_q_rec_des = std::numeric_limits<double>::quiet_NaN();
	m_rec_su_delay = std::numeric_limits<double>::quiet_NaN();
	m_rec_qf_delay = std::numeric_limits<double>::quiet_NaN();
	m_m_dot_max_frac = std::numeric_limits<double>::quiet_NaN();
	m_A_sf = std::numeric_limits<double>::quiet_NaN();

	m_pipe_loss_per_m = std::numeric_limits<double>::quiet_NaN();
	m_pipe_length_add = std::numeric_limits<double>::quiet_NaN();
	m_pipe_length_mult = std::numeric_limits<double>::quiet_NaN();

	m_id_tube = std::numeric_limits<double>::quiet_NaN();
	m_A_tube = std::numeric_limits<double>::quiet_NaN();
	m_n_t = -1;
	m_n_flux_x = 0;
	m_n_flux_y = 0;

	m_T_out_target = std::numeric_limits<double>::quiet_NaN();
	m_hel_stow_deploy = std::numeric_limits<double>::quiet_NaN();

	error_msg = "";
	m_mat_tube = -1;

	m_A_node = std::numeric_limits<double>::quiet_NaN();

	m_Q_dot_piping_loss = std::numeric_limits<double>::quiet_NaN();
	m_m_dot_htf_max = std::numeric_limits<double>::quiet_NaN();

	m_itermode = -1;
	m_od_control = std::numeric_limits<double>::quiet_NaN();
	m_eta_field_iter_prev = std::numeric_limits<double>::quiet_NaN();
	m_tol_od = std::numeric_limits<double>::quiet_NaN();
	m_m_dot_htf_des = std::numeric_limits<double>::quiet_NaN();
	m_q_dot_inc_min = std::numeric_limits<double>::quiet_NaN();

	m_mode = -1;
	m_mode_prev = -1;
	m_E_su = std::numeric_limits<double>::quiet_NaN();
	m_E_su_prev = std::numeric_limits<double>::quiet_NaN();
	m_t_su = std::numeric_limits<double>::quiet_NaN();
	m_t_su_prev = std::numeric_limits<double>::quiet_NaN();

	m_m_mixed = std::numeric_limits<double>::quiet_NaN();
	m_LoverD = std::numeric_limits<double>::quiet_NaN();
	m_RelRough = std::numeric_limits<double>::quiet_NaN();


	m_ncall = -1;
}

void C_gen3gas_receiver::init()
{
	ambient_air.SetFluid(ambient_air.Air);

	// Declare instance of fluid class for FIELD fluid
	field_htfProps.SetFluid(HTFProperties::CO2);
	
	
	// Declare instance of htf class for receiver tube material
	if( m_mat_tube == HTFProperties::Stainless_AISI316 || m_mat_tube == HTFProperties::T91_Steel )
	{
		if( !tube_material.SetFluid(m_mat_tube) )
		{
			throw(C_csp_exception("Tube material code not recognized", "Gas receiver"));
		}
	}
	else if( m_mat_tube == HTFProperties::User_defined )
	{
		throw(C_csp_exception("Receiver material currently does not accept user defined properties", "Gas receiver"));
	}
	else
	{
		error_msg = util::format("Receiver material code, %d, is not recognized", m_mat_tube);
		throw(C_csp_exception(error_msg, "Gas receiver"));
	}

	// Unit Conversions
	m_od_tube /= 1.E3;			//[m] Convert from input in [mm]
	m_th_tube /= 1.E3;			//[m] Convert from input in [mm]
	m_T_hot_des += 273.15;	//[K] Convert from input in [C]
	m_T_cold_des += 273.15;	//[K] Convert from input in [C]
	m_q_rec_des *= 1.E6;		//[W] Convert from input in [MW]

	m_id_tube = m_od_tube - 2 * m_th_tube;			//[m] Inner diameter of receiver tube
	m_A_tube = CSP::pi*m_od_tube / 2.0*m_h_rec;	//[m^2] Outer surface area of each tube
	m_n_t = (int)(m_w_rec / m_od_tube);	// The number of tubes per panel
	
	m_n_lines = 2;		//this design uses 2 collection headers for each centerline distributor

	m_A_node = m_h_rec * m_w_rec / 2. / m_n_flux_x ; //[m^2] The area associated with each node

	m_mode = C_csp_collector_receiver::OFF;					//[-] 0 = requires startup, 1 = starting up, 2 = running
	m_itermode = 1;			//[-] 1: Solve for design temp, 2: solve to match mass flow restriction
	m_od_control = 1.0;			//[-] Additional defocusing for over-design conditions
	m_tol_od = 0.001;		//[-] Tolerance for over-design iteration

	double c_htf_des = field_htfProps.Cp((m_T_hot_des + m_T_cold_des) / 2.0)*1000.0;		//[J/kg-K] Specific heat at design conditions
	m_m_dot_htf_des = m_q_rec_des / (c_htf_des*(m_T_hot_des - m_T_cold_des));					//[kg/s]
	double eta_therm_des = 0.9;
	m_q_dot_inc_min = m_q_rec_des * m_f_rec_min / eta_therm_des;	//[W] Minimum receiver thermal power

	if (m_m_dot_max_frac != m_m_dot_max_frac)
	{
		// if max frac not set, then max mass flow (absolute) needs to be defined
		if (m_m_dot_htf_max != m_m_dot_htf_max)
		{
			throw(C_csp_exception("maximum rec htf mass flow rate not defined", "Gas receiver"));
		}
		m_m_dot_htf_max /= 3600.0;	//[kg/s] Convert from input in [kg/hr]
	}
	m_m_dot_htf_max = m_m_dot_max_frac * m_m_dot_htf_des;	//[kg/s]

	m_mode_prev = m_mode;
	m_E_su_prev = m_q_rec_des * m_rec_qf_delay;	//[W-hr] Startup energy
	m_t_su_prev = m_rec_su_delay;				//[hr] Startup time requirement
	m_eta_field_iter_prev = 1.0;				//[-] Set to largest possible value

	m_T_out_target += 273.15;			//[K] convert from C
	
	// 8.10.2015 twn: Calculate constant thermal losses to the environment
	if(m_pipe_loss_per_m > 0.0 && m_pipe_length_mult > 0.0)
		m_Q_dot_piping_loss = m_pipe_loss_per_m*(m_h_tower*m_pipe_length_mult + m_pipe_length_add);		//[Wt]
	else
		m_Q_dot_piping_loss = 0.0;


	// *******************************************************************

	m_q_dot_inc.resize(m_n_flux_x);
	m_q_dot_inc.fill(0.0);

	m_T_s_guess.resize(m_n_flux_x);
	m_T_s_guess.fill(0.0);
	m_T_s.resize(m_n_flux_x);
	m_T_s.fill(0.0);

	m_T_panel_out_guess.resize(m_n_flux_x);
	m_T_panel_out.resize(m_n_flux_x);
	m_T_panel_out_guess.fill(0.0);
	m_T_panel_out.fill(0.0);

	m_T_panel_ave.resize(m_n_flux_x);
	m_T_panel_ave.fill(0.0);
	m_T_panel_ave_guess.resize(m_n_flux_x);
	m_T_panel_ave_guess.fill(0.0);

	m_T_film.resize(m_n_flux_x);
	m_T_film.fill(0.0);

	m_q_dot_conv.resize(m_n_flux_x);
	m_q_dot_conv.fill(0.0);

	m_q_dot_rad.resize(m_n_flux_x);
	m_q_dot_rad.fill(0.0);

	m_q_dot_loss.resize(m_n_flux_x);
	m_q_dot_loss.fill(0.0);

	m_q_dot_abs.resize(m_n_flux_x);
	m_q_dot_abs.fill(0.0);

	m_m_mixed = 3.2;	//[-] Exponential for calculating mixed convection

	m_LoverD = m_h_rec / m_id_tube;
	m_RelRough = (4.5e-5) / m_id_tube;	//[-] Relative roughness of the tubes. http:www.efunda.com/formulae/fluids/roughness.cfm


	m_ncall = -1;

	return;
}

void C_gen3gas_receiver::call(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	const C_gen3gas_receiver::S_inputs &inputs,
	const C_csp_solver_sim_info &sim_info)
{
	// Increase call-per-timestep counter
	// Converge() sets it to -1, so on first call this line will adjust it = 0
	m_ncall++;
	
	// Get inputs
	double field_eff = inputs.m_field_eff;					//[-]
	const util::matrix_t<double> *flux_map_input = inputs.m_flux_map_input;
		// When this function is called from TCS solver, input_operation_mode should always be == 2
	int input_operation_mode = inputs.m_input_operation_mode;

	if(input_operation_mode < C_csp_collector_receiver::OFF || input_operation_mode > C_csp_collector_receiver::STEADY_STATE)
	{
		error_msg = util::format("Input operation mode must be either [0,1,2], but value is %d", input_operation_mode);
		throw(C_csp_exception(error_msg, "Gas receiver timestep performance call"));
	}

	// Get sim info 
	double step = sim_info.ms_ts.m_step;			//[s]
	double time = sim_info.ms_ts.m_time;	//[s]

	// Get applicable htf state info
	double T_cold_in = htf_state_in.m_temp;		//[C]

	// Complete necessary conversions/calculations of input variables
	T_cold_in += 273.15;				//[K] Cold salt inlet temp, convert from C
	double P_amb = weather.m_pres*100.0;	//[Pa] Ambient pressure, convert from mbar
	double hour = time / 3600.0;			//[hr] Hour of the year
	double T_dp = weather.m_tdew + 273.15;	//[K] Dewpoint temperature, convert from C
	double T_amb = weather.m_tdry + 273.15;	//[K] Dry bulb temperature, convert from C
	// **************************************************************************************

	// Read in remaining weather inputs from weather output structure
	double zenith = weather.m_solzen;
	double azimuth = weather.m_solazi;
	double v_wind_10 = weather.m_wspd;
	double I_bn = weather.m_beam;


	int n_flux_y = (int)flux_map_input->nrows();
	if(n_flux_y > 1)
	{
		error_msg = util::format("The Molten Salt External Receiver (Type222) model does not currently support 2-dimensional "
			"flux maps. The flux profile in the vertical dimension will be averaged. NY=%d", n_flux_y);
		csp_messages.add_message(C_csp_messages::WARNING, error_msg);
	}
	int n_flux_x = (int)flux_map_input->ncols();
	m_flux_in.resize(n_flux_x);

	double T_sky = CSP::skytemp(T_amb, T_dp, hour);

	// Set current timestep stored values to NaN so we know that code solved for them
	m_mode = -1;
	m_E_su = std::numeric_limits<double>::quiet_NaN();
	m_t_su = std::numeric_limits<double>::quiet_NaN();

	m_itermode = 1;

	double v_wind = log((m_h_tower + m_h_rec / 2) / 0.003) / log(10.0 / 0.003)*v_wind_10;

	double c_p_coolant, rho_coolant, f, u_coolant, q_conv_sum, q_rad_sum, q_dot_inc_sum;
	c_p_coolant = rho_coolant = f = u_coolant = q_conv_sum = q_rad_sum = q_dot_inc_sum = std::numeric_limits<double>::quiet_NaN();
	double eta_therm, m_dot_tot, T_hot_guess, m_dot_salt_tot_ss;
	eta_therm = m_dot_tot = T_hot_guess = m_dot_salt_tot_ss = std::numeric_limits<double>::quiet_NaN();
	bool rec_is_off = false;
	bool rec_is_defocusing = false;
	double field_eff_adj = 0.0;

	// ************* Outputs for ISCC model ****************
	double q_thermal_ss = 0.0;
	double f_rec_timestep = 1.0;
	// *****************************************************

	// Do an initial check to make sure the solar position called is valid
	// If it's not, return the output equal to zeros. Also check to make sure
	// the solar flux is at a certain level, otherwise the correlations aren't valid
	if( input_operation_mode == C_csp_collector_receiver::OFF )
	{
		rec_is_off = true;
	}

	if( zenith>(90.0 - m_hel_stow_deploy) || I_bn <= 1.E-6 || (zenith == 0.0 && azimuth == 180.0) )
	{
		m_mode = C_csp_collector_receiver::OFF;
		rec_is_off = true;
	}

	double T_coolant_prop = (m_T_out_target + T_cold_in) / 2.0;		//[K] The temperature at which the coolant properties are evaluated. Validated as constant (mjw)
	c_p_coolant = field_htfProps.Cp(T_coolant_prop)*1000.0;						//[J/kg-K] Specific heat of the coolant

	double m_dot_max = m_m_dot_htf_max;

	double q_abs_sum = 0.0;
	double err_od = 999.0;	// Reset error before iteration

	if (field_eff < m_eta_field_iter_prev && m_od_control < 1.0)
	{	// Suggests controller applied defocus, so reset *controller* defocus
		m_od_control = fmin(m_od_control + (1.0 - field_eff / m_eta_field_iter_prev), 1.0);
	}

	do
	{
		if( rec_is_off )
			break;

		field_eff_adj = field_eff*m_od_control;

		// Get the values of the flux from the fluxmap and store them  as flux_in(col, row)
		if( I_bn > 1.0 )
		{
			double qspec = I_bn * field_eff_adj*m_A_sf / 1000. / (m_h_rec*m_w_rec / (double)n_flux_x);
			for( int j = 0; j<n_flux_x; j++ ){
				m_flux_in.at(j) = 0.;
				
				for( int i = 0; i<n_flux_y; i++ )
					m_flux_in.at(j) += (*flux_map_input)(i, j) * qspec;	//[kW/m^2];

				m_q_dot_inc.at(j) = m_flux_in.at(j) * m_A_node;
			}
		}
		else
		{
			m_flux_in.fill(0.0);
		}

		
		q_dot_inc_sum = 0.0;
		for( int i = 0; i < m_n_flux_x; i++ )
			q_dot_inc_sum += m_q_dot_inc.at(i);		//[kW] Total power absorbed by receiver

		// Set guess values
		m_T_s_guess.fill(m_T_out_target);		//[K] Guess the temperature for the surface nodes
		m_T_panel_out_guess.fill(T_cold_in);	//[K] Guess values for the fluid temp coming out of the control volume

		double c_guess = field_htfProps.Cp((m_T_out_target + T_cold_in) / 2.0);	//[kJ/kg-K] Estimate the specific heat of the fluid in receiver
		double m_dot_guess = std::numeric_limits<double>::quiet_NaN();
		if( I_bn > 1.E-6 )
		{
			double q_guess = 0.5*q_dot_inc_sum;		//[kW] Estimate the thermal power produced by the receiver				
			m_dot_guess = q_guess / (c_guess*(m_T_out_target - T_cold_in)*m_n_lines);	//[kg/s] Mass flow rate for each flow path
		}
		else	// The tower recirculates at night (based on earlier conditions)
		{
			// Enter recirculation mode, where inlet/outlet temps switch
			m_T_out_target = T_cold_in;
			T_cold_in = m_T_s_guess.at(0);		//T_s_guess is set to T_salt_hot before, so this just completes 
			m_dot_guess = -3500.0 / (c_guess*(m_T_out_target - T_cold_in) / 2.0);
		}
		T_hot_guess = 9999.9;		//[K] Initial guess value for error calculation
		double err = -999.9;					//[-] Relative outlet temperature error
		double tol = 0.001;

		//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		//                            ITERATION STARTS HERE
		//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		int qq_max = 50;
		double m_dot = std::numeric_limits<double>::quiet_NaN();
		int qq = 0;
		q_abs_sum = 0.0;

		while( fabs(err) > tol )
		{
			qq++;

			// if the problem fails to converge after 50 iterations, then the power is likely negligible and
			// ..the zero set can be returned
			if( qq > qq_max )
			{
				m_mode = C_csp_collector_receiver::OFF;  // Set the startup mode
				rec_is_off = true;
				break;
			}

			m_dot = m_dot_guess;

			for( int i = 0; i < m_n_flux_x; i++ )
			{
				m_T_s.at(i) = m_T_s_guess.at(i);
				m_T_panel_out.at(i) = m_T_panel_out_guess.at(i);
				// Now do the actual calculations
				m_T_panel_ave.at(i) = (T_cold_in + m_T_panel_out.at(i)) / 2.0;		//[K] The average coolant temperature in each control volume
				m_T_film.at(i) = (m_T_s.at(i) + T_amb) / 2.0;					//[K] Film temperature
			}

			// Calculate the average surface temperature
			double T_s_sum = 0.0;
			for( int i = 0; i < m_n_flux_x; i++ )
				T_s_sum += m_T_s.at(i);
			double T_film_ave = (T_amb + m_T_out_target) / 2.0;

			// Convective coefficient for external forced convection using Siebers & Kraabel
			double k_film = ambient_air.cond(T_film_ave);				//[W/m-K] The conductivity of the ambient air
			double mu_film = ambient_air.visc(T_film_ave);			//[kg/m-s] Dynamic viscosity of the ambient air
			double rho_film = ambient_air.dens(T_film_ave, P_amb);	//[kg/m^3] Density of the ambient air
			double c_p_film = ambient_air.Cp(T_film_ave);				//[kJ/kg-K] Specific heat of the ambient air
			double Re_for = rho_film*v_wind*m_w_rec / mu_film;			//[-] Reynolds number
			double ksD = (m_od_tube / 2.0) / m_w_rec;						//[-] The effective roughness of the cylinder [Siebers, Kraabel 1984]
			double Nusselt_for = CSP::Nusselt_FC(ksD, Re_for);		//[-] S&K
			double h_for = Nusselt_for*k_film / m_w_rec*m_hl_ffact;		//[W/m^2-K] Forced convection heat transfer coefficient

			// Convection coefficient for external natural convection using Siebers & Kraabel
			// Note: This relationship applies when the surrounding properties are evaluated at ambient conditions [S&K]
			double beta = 1.0 / T_amb;												//[1/K] Volumetric expansion coefficient
			double nu_amb = ambient_air.visc(T_amb) / ambient_air.dens(T_amb, P_amb);	//[m^2/s] Kinematic viscosity		


			for( int i = 0; i < m_n_flux_x; i++ )
			{
				// Natural convection
				double Gr_nat = fmax(0.0, CSP::grav*beta*(m_T_s.at(i) - T_amb)*pow(m_h_rec, 3) / pow(nu_amb, 2));	//[-] Grashof Number at ambient conditions
				double Nusselt_nat = 0.098*pow(Gr_nat, (1.0 / 3.0))*pow(m_T_s.at(i) / T_amb, -0.14);					//[-] Nusselt number
				double h_nat = Nusselt_nat*ambient_air.cond(T_amb) / m_h_rec*m_hl_ffact;							//[W/m^-K] Natural convection coefficient
				// Mixed convection
				double h_mixed = pow((pow(h_for, m_m_mixed) + pow(h_nat, m_m_mixed)), 1.0 / m_m_mixed)*4.0;			//(4.0) is a correction factor to match convection losses at Solar II (correspondance with G. Kolb, SNL)
				m_q_dot_conv.at(i) = h_mixed*m_A_node*(m_T_s.at(i) - m_T_film.at(i));							//[W] Convection losses per node
				// Radiation from the receiver - Calculate the radiation node by node
				m_q_dot_rad.at(i) = 0.5*CSP::sigma*m_epsilon*m_A_node*(2.0*pow(m_T_s.at(i), 4) - pow(T_amb, 4) - pow(T_sky, 4))*m_hl_ffact;	//[W] Total radiation losses per node
				m_q_dot_loss.at(i) = m_q_dot_rad.at(i) + m_q_dot_conv.at(i);			//[W] Total overall losses per node
				m_q_dot_abs.at(i) = m_q_dot_inc.at(i)*1000.0 - m_q_dot_loss.at(i);	//[W] Absorbed flux at each node
				// Calculate the temperature drop across the receiver tube wall... assume a cylindrical thermal resistance
				double T_wall = (m_T_s.at(i) + m_T_panel_ave.at(i)) / 2.0;				//[K] The temperature at which the conductivity of the wall is evaluated
				double k_tube = tube_material.cond(T_wall);								//[W/m-K] The conductivity of the wall
				double R_tube_wall = m_th_tube / (k_tube*m_A_node);	//[K/W] The thermal resistance of the wall
				// Calculations for the inside of the tube						
				double mu_coolant = field_htfProps.visc(T_coolant_prop);					//[kg/m-s] Absolute viscosity of the coolant
				double k_coolant = field_htfProps.cond(T_coolant_prop);					//[W/m-K] Conductivity of the coolant
				rho_coolant = field_htfProps.dens(T_coolant_prop, 1.0);			//[kg/m^3] Density of the coolant

				u_coolant = m_dot / (m_n_t*rho_coolant*pow((m_id_tube / 2.0), 2)*CSP::pi);	//[m/s] Average velocity of the coolant through the receiver tubes
				double Re_inner = rho_coolant*u_coolant*m_id_tube / mu_coolant;				//[-] Reynolds number of internal flow
				double Pr_inner = c_p_coolant*mu_coolant / k_coolant;							//[-] Prandtl number of internal flow
				double Nusselt_t;
				CSP::PipeFlow(Re_inner, Pr_inner, m_LoverD, m_RelRough, Nusselt_t, f);
				if( Nusselt_t <= 0.0 )
				{
					m_mode = C_csp_collector_receiver::OFF;		// Set the startup mode
					rec_is_off = true;
					break;
				}
				double h_inner = Nusselt_t*k_coolant / m_id_tube;								//[W/m^2-K] Convective coefficient between the inner tube wall and the coolant
				double R_conv_inner = 1.0 / (h_inner*CSP::pi*m_id_tube / 2.0*m_h_rec*m_n_t);	//[K/W] Thermal resistance associated with this value


				m_T_panel_out_guess.at(i) = T_cold_in + m_q_dot_abs.at(i) / (m_dot*c_p_coolant);	//[K] Energy balance for each node																																																
				m_T_panel_ave_guess.at(i) = (m_T_panel_out_guess.at(i) + T_cold_in) / 2.0;				//[K] Panel average temperature
				m_T_s_guess.at(i) = m_T_panel_ave_guess.at(i) + m_q_dot_abs.at(i)*(R_conv_inner + R_tube_wall);			//[K] Surface temperature based on the absorbed heat

				if( m_T_s_guess.at(i) < 1.0 )
				{
					m_mode = C_csp_collector_receiver::OFF;  // Set the startup mode
					rec_is_off = true;
				}
			}

			if( rec_is_off )
				break;

			q_conv_sum = 0.0; q_rad_sum = 0.0; //q_inc_sum = 0.0;
			q_abs_sum = 0.0;
			for( int i = 0; i < m_n_flux_x; i++ )
			{
				q_conv_sum += m_q_dot_conv.at(i);
				q_rad_sum += m_q_dot_rad.at(i);
				//q_inc_sum += m_q_dot_inc.at(i,0);
				q_abs_sum += m_q_dot_abs.at(i);
			}

			T_hot_guess = 0.0;
			for( int j = 0; j < m_n_flux_x; j++ )
				T_hot_guess += m_T_panel_out_guess.at(j);		//[K] Update the calculated hot salt outlet temp
			T_hot_guess /= (double)m_n_flux_x;

			
			//Calculate outlet temperature after piping losses
			if( m_Q_dot_piping_loss > 0.0 )
			{
				double tmp = m_dot*m_n_lines;		//[kg/s]
				double delta_T_piping = m_Q_dot_piping_loss / (tmp*c_p_coolant);	//[K]
				T_hot_guess = T_hot_guess - m_Q_dot_piping_loss/(tmp*c_p_coolant);	//[K]
			}

			//calculate convergence error
			err = (T_hot_guess - m_T_out_target) / m_T_out_target;

			if( fabs(err) > tol )
			{
				m_dot_guess = (q_abs_sum - m_Q_dot_piping_loss) / (m_n_lines*c_p_coolant*(m_T_out_target - T_cold_in));			//[kg/s]

				if( m_dot_guess < 1.E-5 )
				{
					m_mode = C_csp_collector_receiver::OFF;				//[-] Set the startup mode
					rec_is_off = true;
				}
			}
		}

		if( rec_is_off )
			break;

		// Now we can calculate some of the parasitics associated with pumping the coolant fluid
		// Calculating the pressure drop across the receiver panels
		m_dot_tot = m_dot*m_n_lines;
		double m_dot_tube = m_dot / (double)m_n_t;		//[kg/s] The mass flow through each individual tube

		// Limit the HTF mass flow rate to the maximum, if needed
		if( (m_dot_tot > m_dot_max) || m_itermode == 2 )
		{
			err_od = (m_dot_tot - m_dot_max) / m_dot_max;
			if( err_od < m_tol_od )
			{
				m_itermode = 1;
				m_od_control = 1.0;
				rec_is_defocusing = false;
			}
			else
			{
				m_od_control = m_od_control*pow((m_dot_max / m_dot_tot), 0.8);	//[-] Adjust the over-design defocus control by modifying the current value
				m_itermode = 2;
				rec_is_defocusing = true;
				// GOTO 15
			}
		}
	} while( rec_is_defocusing );

	// 8.12.2015 twn: not using eta_therm in iteration loop - so move this calculation after loop
	if( q_dot_inc_sum > 0.0 )
		eta_therm = q_abs_sum / (q_dot_inc_sum*1000.0);
	else
		eta_therm = 0.0;

	double DELTAP, Pres_D, W_dot_pump, q_thermal, q_startup;
	DELTAP = Pres_D = W_dot_pump = q_thermal = q_startup = std::numeric_limits<double>::quiet_NaN();

	q_startup = 0.0;

	double time_required_su = step/3600.0;

	if( !rec_is_off )
	{
		m_dot_salt_tot_ss = m_dot_tot;

		switch( input_operation_mode )
		{
		case C_csp_collector_receiver::STARTUP:
			{
				double time_require_su_energy = m_E_su_prev / (m_dot_tot*c_p_coolant*(T_hot_guess - T_cold_in));	//[hr]
				double time_require_su_ramping = m_t_su_prev;

				double time_required_max = fmax(time_require_su_energy, time_require_su_ramping);	//[hr]

				double time_step_hrs = step / 3600.0;		//[hr]

				if( time_required_max  > time_step_hrs )		// Can't completely startup receiver in maximum allowable timestep
				{											// Need to advance timestep and try again
					time_required_su = time_step_hrs;		
					m_mode = C_csp_collector_receiver::STARTUP;
					q_startup = m_dot_tot*c_p_coolant*(T_hot_guess - T_cold_in)*step / 3600.0;
				}
				else
				{
					time_required_su = time_required_max;		//[hr]
					m_mode = C_csp_collector_receiver::ON;

					double q_startup_energy_req = m_E_su_prev;	//[W-hr]
					double q_startup_ramping_req = m_dot_tot*c_p_coolant*(T_hot_guess - T_cold_in)*m_t_su;	//[W-hr]
					q_startup = fmax(q_startup_energy_req, q_startup_ramping_req);
				}

				m_E_su = fmax(0.0, m_E_su_prev - m_dot_tot*c_p_coolant*(T_hot_guess - T_cold_in)*step / 3600.0);
				m_t_su = fmax(0.0, m_t_su_prev - step / 3600.0);
			}

			rec_is_off = true;

			break;

		case C_csp_collector_receiver::ON:
			
			if( m_E_su_prev > 0.0 || m_t_su_prev > 0.0 )
			{
				
				m_E_su = fmax(0.0, m_E_su_prev - m_dot_tot*c_p_coolant*(T_hot_guess - T_cold_in)*step / 3600.0);	//[W-hr]
				m_t_su = fmax(0.0, m_t_su_prev - step / 3600.0);	//[hr]

				if( m_E_su + m_t_su > 0.0 )
				{
					m_mode = C_csp_collector_receiver::STARTUP;		// If either are greater than 0, we're staring up but not finished
					
					// 4.28.15 twn: Startup energy also needs to consider energy consumed during time requirement, if that is greater than energy requirement
						//q_startup = (m_E_su_prev - m_E_su) / (step / 3600.0)*1.E-6;
					q_startup = m_dot_tot*c_p_coolant*(T_hot_guess - T_cold_in)*step / 3600.0;

					rec_is_off = true;
					f_rec_timestep = 0.0;
				}
				else
				{
					m_mode = C_csp_collector_receiver::ON;

					double q_startup_energy_req = m_E_su_prev;	//[W-hr]
					double q_startup_ramping_req = m_dot_tot*c_p_coolant*(T_hot_guess - T_cold_in)*m_t_su;	//[W-hr]
					q_startup = fmax(q_startup_energy_req, q_startup_ramping_req);

					// Adjust the available mass flow to reflect startup
					m_dot_tot = fmin((1.0 - m_t_su_prev / (step / 3600.0))*m_dot_tot, m_dot_tot - m_E_su_prev / ((step / 3600.0)*c_p_coolant*(T_hot_guess - T_cold_in)));
					f_rec_timestep = fmax(0.0, fmin(1.0 - m_t_su_prev / (step / 3600.0), 1.0 - m_E_su_prev / (m_dot_tot*c_p_coolant*(T_hot_guess - T_cold_in))));
				}
					//4.28.15 twn: Startup energy needs to consider
				//q_startup = (m_E_su_prev - m_E_su) / (step / 3600.0)*1.E-6;
			}
			else
			{
				m_E_su = m_E_su_prev;
				m_t_su = m_t_su_prev;
				m_mode = C_csp_collector_receiver::ON;
				q_startup = 0.0;

				q_thermal = m_dot_tot*c_p_coolant*(T_hot_guess - T_cold_in);

				//if( q_thermal < m_q_rec_min )
				if(q_dot_inc_sum*1.E3 < m_q_dot_inc_min)
				{
					// If output here is less than specified allowed minimum, then need to shut off receiver
					m_mode = C_csp_collector_receiver::OFF;

					// Include here outputs that are ONLY set to zero if receiver completely off, and not attempting to start-up
					W_dot_pump = 0.0;
					// Pressure drops
					DELTAP = 0.0; Pres_D = 0.0; u_coolant = 0.0;
				}
			}
			break;

		case C_csp_collector_receiver::STEADY_STATE:

			m_mode = C_csp_collector_receiver::STEADY_STATE;
			f_rec_timestep = 1.0;

			break;
		
		}	// End switch() on input_operation_mode

		// Pressure drop calculations
        calc_receiver_presdrop(rho_coolant, m_dot_tot, f, Pres_D);

		q_thermal = m_dot_tot*c_p_coolant*(T_hot_guess - T_cold_in);
		q_thermal_ss = m_dot_salt_tot_ss*c_p_coolant*(T_hot_guess - T_cold_in);

		// After convergence, determine whether the mass flow rate falls below the lower limit
		if(q_dot_inc_sum*1.E3 < m_q_dot_inc_min)
		{
			// GOTO 900
			// Steady State always reports q_thermal (even when much less than min) because model is letting receiver begin startup with this energy
			// Should be a way to communicate to controller that q_thermal is less than q_min without losing this functionality
			if(m_mode != C_csp_collector_receiver::STEADY_STATE || m_mode_prev == C_csp_collector_receiver::ON)
				rec_is_off = true;
		}
	}
	else
	{	// If receiver was off BEFORE startup deductions
		m_mode = C_csp_collector_receiver::OFF;

		// Include here outputs that are ONLY set to zero if receiver completely off, and not attempting to start-up
		W_dot_pump = 0.0;
		// Pressure drops
		DELTAP = 0.0; Pres_D = 0.0; u_coolant = 0.0;
	}

	if( rec_is_off )
	{
		// 900 continue	// Receiver isn't producing usable energy
		m_dot_tot = 0.0; eta_therm = 0.0; /*W_dot_pump = 0.0;*/
		q_conv_sum = 0.0; q_rad_sum = 0.0; m_T_s.fill(0.0); q_thermal = 0.0;
		// Set the receiver outlet temperature equal to the inlet design temperature
		T_hot_guess = m_T_cold_des;
		q_dot_inc_sum = 0.0;
		// Pressure drops
		/*DELTAP = 0.0; Pres_D = 0.0; u_coolant = 0.0;*/
		// Set receiver startup energy to 0
		// q_startup = 0.0;
		// ISCC outputs
		m_dot_salt_tot_ss = 0.0; f_rec_timestep = 0.0; q_thermal_ss = 0.0;

		// Reset m_od_control
		m_od_control = 1.0;		//[-]
	}

	ms_outputs.m_m_dot_salt_tot = m_dot_tot*3600.0;		//[kg/hr] convert from kg/s
	ms_outputs.m_eta_therm = eta_therm;							//[-] RECEIVER thermal efficiency (includes radiation and convective losses. reflection losses are contained in receiver flux model)
	ms_outputs.m_W_dot_pump = W_dot_pump / 1.E6;				//[MW] convert from W
	ms_outputs.m_q_conv_sum = q_conv_sum / 1.E6;				//[MW] convert from W
	ms_outputs.m_q_rad_sum = q_rad_sum / 1.E6;					//[MW] convert from W
	ms_outputs.m_Q_thermal = q_thermal / 1.E6;					//[MW] convert from W
	ms_outputs.m_T_salt_hot = T_hot_guess - 273.15;		//[C] convert from K
	ms_outputs.m_field_eff_adj = field_eff_adj;					//[-]
	ms_outputs.m_component_defocus = m_od_control;				//[-]
	ms_outputs.m_q_dot_rec_inc = q_dot_inc_sum / 1.E3;			//[MW] convert from kW
	ms_outputs.m_q_startup = q_startup/1.E6;					//[MW-hr] convert from W-hr
	ms_outputs.m_dP_receiver = DELTAP / 1.E5;	//[bar] receiver pressure drop, convert from Pa
	ms_outputs.m_dP_total = Pres_D*10.0;						//[bar] total pressure drop, convert from MPa
	ms_outputs.m_vel_htf = u_coolant;							//[m/s]
	ms_outputs.m_T_salt_cold = T_cold_in - 273.15;			//[C] convert from K
	ms_outputs.m_m_dot_ss = m_dot_salt_tot_ss*3600.0;			//[kg/hr] convert from kg/s
	ms_outputs.m_q_dot_ss = q_thermal_ss / 1.E6;				//[MW] convert from W
	ms_outputs.m_f_timestep = f_rec_timestep;					//[-]
	ms_outputs.m_time_required_su = time_required_su*3600.0;	//[s], convert from hr in code
	if(q_thermal > 0.0)
		ms_outputs.m_q_dot_piping_loss = m_Q_dot_piping_loss/1.E6;	//[MWt]
	else
		ms_outputs.m_q_dot_piping_loss = 0.0;		//[MWt]

	m_eta_field_iter_prev = field_eff;	//[-]
}

void C_gen3gas_receiver::off(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	const C_csp_solver_sim_info &sim_info)
{
	// Don't currently need *any* of these inputs, but if we add recirculation or thermal capacitance it would be helpful to have in place
	m_mode = C_csp_collector_receiver::OFF;

	// Assuming no night recirculation, so... these should be zero
	ms_outputs.m_m_dot_salt_tot = 0.0;		//[kg/hr] convert from kg/s
	ms_outputs.m_eta_therm = 0.0;			//[-] RECEIVER thermal efficiency (includes radiation and convective losses. reflection losses are contained in receiver flux model)
	ms_outputs.m_W_dot_pump = 0.0;			//[MW] convert from W
	ms_outputs.m_q_conv_sum = 0.0;			//[MW] convert from W
	ms_outputs.m_q_rad_sum = 0.0;			//[MW] convert from W
	ms_outputs.m_Q_thermal = 0.0;			//[MW] convert from W
	ms_outputs.m_T_salt_hot = 0.0;			//[C] convert from K
	ms_outputs.m_field_eff_adj = 0.0;		//[-]
	ms_outputs.m_component_defocus = 1.0;	//[-]
	ms_outputs.m_q_dot_rec_inc = 0.0;		//[MW] convert from kW
	ms_outputs.m_q_startup = 0.0;			//[MW-hr] convert from W-hr
	ms_outputs.m_dP_receiver = 0.0;			//[bar] receiver pressure drop, convert from Pa
	ms_outputs.m_dP_total = 0.0;			//[bar] total pressure drop, convert from MPa
	ms_outputs.m_vel_htf = 0.0;				//[m/s]
	ms_outputs.m_T_salt_cold = 0.0;			//[C] convert from K
	ms_outputs.m_m_dot_ss = 0.0;			//[kg/hr] convert from kg/s
	ms_outputs.m_q_dot_ss = 0.0;			//[MW] convert from W
	ms_outputs.m_f_timestep = 0.0;			//[-]
	ms_outputs.m_time_required_su = sim_info.ms_ts.m_step;	//[s], convert from hr in code
	ms_outputs.m_q_dot_piping_loss = 0.0;	//[MWt]
	
	return;
}

void C_gen3gas_receiver::converged()
{
	// Check HTF props?
	//!MJW 9.8.2010 :: Call the property range check subroutine with the inlet and outlet HTF temps to make sure they're in the valid range
	//call check_htf(Coolant,T_salt_hot)
	//call check_htf(Coolant,T_salt_cold)

	if( m_mode == C_csp_collector_receiver::STEADY_STATE )
	{
		throw(C_csp_exception("Receiver should only be run at STEADY STATE mode for estimating output. It must be run at a different mode before exiting a timestep",
			"Gas receiver converged method"));
	}

	if( m_mode == C_csp_collector_receiver::OFF )
	{
		m_E_su = m_q_rec_des * m_rec_qf_delay;
		m_t_su = m_rec_su_delay;
	}

	m_mode_prev = m_mode;
	m_E_su_prev = m_E_su;
	m_t_su_prev = m_t_su;

	m_itermode = 1;
	m_od_control = 1.0;
	m_eta_field_iter_prev = 1.0;		//[-]

	m_ncall = -1;
}

int C_gen3gas_receiver::get_operating_state()
{
	return m_mode_prev;
}

void C_gen3gas_receiver::clear_outputs()
{
	ms_outputs.m_m_dot_salt_tot = 
		ms_outputs.m_eta_therm = 
		ms_outputs.m_W_dot_pump = 
		ms_outputs.m_q_conv_sum = 
		ms_outputs.m_q_rad_sum = 
		ms_outputs.m_Q_thermal =
		ms_outputs.m_T_salt_hot = 
		ms_outputs.m_field_eff_adj = 
		ms_outputs.m_component_defocus =
		ms_outputs.m_q_dot_rec_inc = 
		ms_outputs.m_q_startup = 
		ms_outputs.m_dP_receiver = 
		ms_outputs.m_dP_total =
		ms_outputs.m_vel_htf = 
		ms_outputs.m_T_salt_cold = 
		ms_outputs.m_m_dot_ss = 
		ms_outputs.m_q_dot_ss = 
		ms_outputs.m_f_timestep = std::numeric_limits<double>::quiet_NaN();
}

void C_gen3gas_receiver::calc_receiver_presdrop(double rho_f, double mdot, double ffact, double &PresDrop_calc)
{

    // Pressure drop calculations

}

HTFProperties *C_gen3gas_receiver::get_htf_property_object()
{
    return &field_htfProps;
}