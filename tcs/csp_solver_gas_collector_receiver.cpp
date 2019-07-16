/*******************************************************************************************************
*  Copyright 2017 Alliance for Sustainable Energy, LLC
*
*  NOTICE: This software was developed at least in part by Alliance for Sustainable Energy, LLC
*  (�Alliance�) under Contract No. DE-AC36-08GO28308 with the U.S. Department of Energy and the U.S.
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
*  the underlying software originally provided by Alliance as �System Advisor Model� or �SAM�. Except
*  to comply with the foregoing, the terms �System Advisor Model�, �SAM�, or any confusingly similar
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

#include "csp_solver_gas_collector_receiver.h"
#include "sam_csp_util.h"
#include <algorithm>

static C_csp_reported_outputs::S_output_info S_output_info[] =
{
	{C_csp_gas_collector_receiver::E_FIELD_Q_DOT_INC, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_gas_collector_receiver::E_FIELD_ETA_OPT, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_gas_collector_receiver::E_FIELD_ADJUST, C_csp_reported_outputs::TS_WEIGHTED_AVE},

	{C_csp_gas_collector_receiver::E_Q_DOT_INC, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_gas_collector_receiver::E_ETA_THERMAL, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_gas_collector_receiver::E_Q_DOT_THERMAL, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_gas_collector_receiver::E_M_DOT_HTF, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_gas_collector_receiver::E_Q_DOT_STARTUP, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_gas_collector_receiver::E_T_HTF_IN, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_gas_collector_receiver::E_T_HTF_OUT, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_gas_collector_receiver::E_Q_DOT_PIPE_LOSS, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_gas_collector_receiver::E_Q_DOT_LOSS, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	
	csp_info_invalid	
};

C_csp_gas_collector_receiver::C_csp_gas_collector_receiver(
	std::vector<C_pt_sf_perf_interp> & pt_heliostatfields,
	std::vector<C_gen3gas_receiver> & gen3gas_receivers
):
	mc_pt_heliostatfields(pt_heliostatfields),
	mc_gen3gas_receivers(gen3gas_receivers)
{
	mc_reported_outputs.construct(S_output_info);
}

C_csp_gas_collector_receiver::~C_csp_gas_collector_receiver()
{}

void C_csp_gas_collector_receiver::init(const C_csp_collector_receiver::S_csp_cr_init_inputs init_inputs, 
				C_csp_collector_receiver::S_csp_cr_solved_params & solved_params)
{
	
	if(mc_pt_heliostatfields.size() != mc_gen3gas_receivers.size())
		throw(C_csp_exception("The number of heliostat fields and receivers supplied to the model must be equal.", "Gas receiver"));

	
	for (size_t i = 0; i < mc_gen3gas_receivers.size(); i++)
	{
		mc_pt_heliostatfields.at(i).init();
		mc_gen3gas_receivers.at(i).init();
	}

	solved_params.m_T_htf_cold_des = mc_gen3gas_receivers.front().m_T_cold_des;			//[K]
	solved_params.m_q_dot_rec_des = mc_gen3gas_receivers.front().m_q_rec_des / 1.E6;		//[MW]

	solved_params.m_A_aper_total = 0.;
	for (size_t i = 0; i < mc_pt_heliostatfields.size(); i++)
		solved_params.m_A_aper_total += mc_pt_heliostatfields.at(i).ms_params.m_A_sf;			//[m^2]

	return;
}

int C_csp_gas_collector_receiver::get_operating_state()
{
	return mc_gen3gas_receivers.front().get_operating_state();
}

double C_csp_gas_collector_receiver::get_startup_time()
{
    return mc_gen3gas_receivers.front().m_rec_su_delay * 3600.; //sec   
}

double C_csp_gas_collector_receiver::get_startup_energy() //MWh
{
    return mc_gen3gas_receivers.front().m_rec_qf_delay * mc_gen3gas_receivers.front().m_q_rec_des*1.e-6;
}

double C_csp_gas_collector_receiver::get_pumping_parasitic_coef()  //MWe/MWt
{
	return 0.;
}

double C_csp_gas_collector_receiver::get_min_power_delivery()    //MWt
{
    return mc_gen3gas_receivers.front().m_f_rec_min * mc_gen3gas_receivers.front().m_q_rec_des*1.e-6;
}


double C_csp_gas_collector_receiver::get_tracking_power()
{
	double track_tot = 0.;
	for (size_t i = 0; i < mc_pt_heliostatfields.size(); i++)
		track_tot += mc_pt_heliostatfields.at(i).ms_params.m_p_track * mc_pt_heliostatfields.at(i).ms_params.m_N_hel*1.e-3;	//MWe
	return track_tot;
}

double C_csp_gas_collector_receiver::get_col_startup_power()
{
	double start_tot = 0.;
	for (size_t i = 0; i < mc_pt_heliostatfields.size(); i++)
		start_tot += mc_pt_heliostatfields.at(i).ms_params.m_p_start * mc_pt_heliostatfields.at(i).ms_params.m_N_hel*1.e-3;	//MWe
	return start_tot;
}


void C_csp_gas_collector_receiver::call(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	const C_csp_collector_receiver::S_csp_cr_inputs &inputs,
	C_csp_collector_receiver::S_csp_cr_out_solver &cr_out_solver,
	//C_csp_collector_receiver::S_csp_cr_out_report &cr_out_report,
	const C_csp_solver_sim_info &sim_info)
{
	// First call heliostat field class: 'csp_solver_pt_heliostat'
	// Then use its outputs as inputs to receiver class: 'csp_solver_mspt_receiver_222'

	// Set heliostat field call() parameters and solve
	double heliostat_field_control = inputs.m_field_control;
	for (size_t i = 0; i < mc_pt_heliostatfields.size(); i++)
		mc_pt_heliostatfields.at(i).call(weather, heliostat_field_control, sim_info);

	// Get heliostat field outputs and set corresponding receiver inputs
	for (size_t i = 0; i < mc_gen3gas_receivers.size(); i++)
	{
		C_gen3gas_receiver::S_inputs receiver_inputs;
		receiver_inputs.m_field_eff = mc_pt_heliostatfields.at(i).ms_outputs.m_eta_field;
		receiver_inputs.m_input_operation_mode = inputs.m_input_operation_mode;
		receiver_inputs.m_flux_map_input = &mc_pt_heliostatfields.at(i).ms_outputs.m_flux_map_out;
		mc_gen3gas_receivers.at(i).call(weather, htf_state_in, receiver_inputs, sim_info);
	}
		
	// Set collector/receiver parent class outputs and return
	cr_out_solver.m_q_thermal = 0.;				//[MW]
	cr_out_solver.m_q_startup = 0.;				//[MWt-hr]
	cr_out_solver.m_m_dot_salt_tot = 0.;		//[kg/hr]
	cr_out_solver.m_T_salt_hot = 0.;				//[C]
	cr_out_solver.m_component_defocus = 0.;	//[-]
	cr_out_solver.m_W_dot_htf_pump = 0.;			//[MWe]
	cr_out_solver.m_W_dot_col_tracking = 0.;		//[MWe]
	cr_out_solver.m_time_required_su = 0.;	//[s]

	double norm = 1. / (double)mc_gen3gas_receivers.size();

	double q_dot_inc = 0, eta_thermal = 0., q_dot_thermal = 0., m_dot_htf = 0., 
			q_dot_startup = 0., T_htf_in = 0., T_htf_out = 0., q_dot_pipe_loss = 0., 
			q_dot_loss = 0., q_dot_field_inc = 0., eta_field=0., sf_adjust_out=0.;

	for (size_t i = 0; i < mc_gen3gas_receivers.size(); i++)
	{
		cr_out_solver.m_q_thermal += mc_gen3gas_receivers.at(i).ms_outputs.m_Q_thermal;				//[MW]
		cr_out_solver.m_q_startup += mc_gen3gas_receivers.at(i).ms_outputs.m_q_startup;				//[MWt-hr]
		cr_out_solver.m_T_salt_hot += mc_gen3gas_receivers.at(i).ms_outputs.m_T_salt_hot*norm;				//[C]
		cr_out_solver.m_component_defocus += mc_gen3gas_receivers.at(i).ms_outputs.m_component_defocus;	//[-]
		cr_out_solver.m_W_dot_htf_pump += mc_gen3gas_receivers.at(i).ms_outputs.m_W_dot_pump;			//[MWe]
		cr_out_solver.m_W_dot_col_tracking += mc_pt_heliostatfields.at(i).ms_outputs.m_pparasi;		//[MWe]
		cr_out_solver.m_time_required_su = std::fmax(cr_out_solver.m_time_required_su, mc_gen3gas_receivers.at(i).ms_outputs.m_time_required_su);	//[s]

		q_dot_inc += mc_gen3gas_receivers.at(i).ms_outputs.m_q_dot_rec_inc;
		eta_thermal += mc_gen3gas_receivers.at(i).ms_outputs.m_eta_therm*mc_gen3gas_receivers.at(i).ms_outputs.m_Q_thermal;
		q_dot_thermal += mc_gen3gas_receivers.at(i).ms_outputs.m_Q_thermal;
		q_dot_startup += mc_gen3gas_receivers.at(i).ms_outputs.m_q_startup / (mc_gen3gas_receivers.at(i).ms_outputs.m_time_required_su / 3600.0);
		T_htf_out += mc_gen3gas_receivers.at(i).ms_outputs.m_T_salt_hot*norm;
		q_dot_pipe_loss += mc_gen3gas_receivers.at(i).ms_outputs.m_q_dot_piping_loss;
		q_dot_loss += mc_gen3gas_receivers.at(i).ms_outputs.m_q_rad_sum + mc_gen3gas_receivers.at(i).ms_outputs.m_q_conv_sum;
		q_dot_field_inc += mc_pt_heliostatfields.at(i).ms_outputs.m_q_dot_field_inc;
		eta_field += mc_pt_heliostatfields.at(i).ms_outputs.m_eta_field * mc_pt_heliostatfields.at(i).ms_outputs.m_q_dot_field_inc;
		sf_adjust_out += mc_pt_heliostatfields.at(i).ms_outputs.m_sf_adjust_out * mc_pt_heliostatfields.at(i).ms_outputs.m_q_dot_field_inc;
	}
	//The "total" mass flow is the max flow through the system, and this is governed by the first (north) receiver
	m_dot_htf = cr_out_solver.m_m_dot_salt_tot = mc_gen3gas_receivers.front().ms_outputs.m_m_dot_tot;		//[kg/hr]

	eta_thermal /= q_dot_thermal; //weighted efficiency
	eta_field /= q_dot_field_inc;
	sf_adjust_out /= q_dot_field_inc;

	mc_reported_outputs.value(E_FIELD_Q_DOT_INC, q_dot_field_inc);	//[MWt]
	mc_reported_outputs.value(E_FIELD_ETA_OPT, eta_field);			//[-]
	mc_reported_outputs.value(E_FIELD_ADJUST, sf_adjust_out);			//[-]

	mc_reported_outputs.value(E_Q_DOT_INC, q_dot_inc);	//[MWt]
	mc_reported_outputs.value(E_ETA_THERMAL, eta_thermal);		//[-]
	mc_reported_outputs.value(E_Q_DOT_THERMAL, q_dot_thermal);	//[MWt]
	mc_reported_outputs.value(E_M_DOT_HTF, m_dot_htf);	//[kg/hr]
		// If startup, then timestep may have changed (why not report this from 222 in MWt?)
	mc_reported_outputs.value(E_Q_DOT_STARTUP, q_dot_startup);		//[MWt])
	mc_reported_outputs.value(E_T_HTF_IN, htf_state_in.m_temp);									//[C]
	mc_reported_outputs.value(E_T_HTF_OUT, T_htf_out);		//[C]
	mc_reported_outputs.value(E_Q_DOT_PIPE_LOSS, q_dot_pipe_loss);	//[MWt]
    mc_reported_outputs.value(E_Q_DOT_LOSS, q_dot_loss); //MWt
}

void C_csp_gas_collector_receiver::off(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	C_csp_collector_receiver::S_csp_cr_out_solver &cr_out_solver,
	//C_csp_collector_receiver::S_csp_cr_out_report &cr_out_report,
	const C_csp_solver_sim_info &sim_info)
{
	// First call heliostat field class
	// In OFF call, looking specifically for weather STOW parasitics apply
	for (size_t i = 0; i < mc_gen3gas_receivers.size(); i++)
		mc_pt_heliostatfields.at(i).off(sim_info);

	// Set collector/receiver parent class outputs from field model
	cr_out_solver.m_component_defocus = 1.0;	//[-]

	cr_out_solver.m_W_dot_col_tracking = 0.;	//[MWe]
	cr_out_solver.m_q_thermal = 0.;				//[MW]
	cr_out_solver.m_q_startup = 0.;				//[MWt-hr]
	cr_out_solver.m_m_dot_salt_tot = 0.;		//[kg/hr]
	cr_out_solver.m_T_salt_hot = 0.;				//[C]
	cr_out_solver.m_component_defocus = 0.;	//[-]
	cr_out_solver.m_W_dot_htf_pump = 0.;			//[MWe]
	cr_out_solver.m_W_dot_col_tracking = 0.;		//[MWe]
	cr_out_solver.m_time_required_su = 0.;	//[s]

	double norm = 1. / (double)mc_gen3gas_receivers.size();

	double q_dot_inc = 0, eta_thermal = 0., q_dot_thermal = 0., m_dot_htf = 0.,
		q_dot_startup = 0., T_htf_in = 0., T_htf_out = 0., q_dot_pipe_loss = 0.,
		q_dot_loss = 0., q_dot_field_inc = 0., eta_field = 0., sf_adjust_out = 0.;

	// Set collector/receiver parent class outputs from field model
	for (size_t i = 0; i < mc_gen3gas_receivers.size(); i++)
	{
		mc_gen3gas_receivers.at(i).off(weather, htf_state_in, sim_info);

		cr_out_solver.m_q_thermal += mc_gen3gas_receivers.at(i).ms_outputs.m_Q_thermal;				//[MW]
		cr_out_solver.m_q_startup += mc_gen3gas_receivers.at(i).ms_outputs.m_q_startup;				//[MWt-hr]
		cr_out_solver.m_T_salt_hot += mc_gen3gas_receivers.at(i).ms_outputs.m_T_salt_hot*norm;				//[C]
		cr_out_solver.m_W_dot_htf_pump += mc_gen3gas_receivers.at(i).ms_outputs.m_W_dot_pump;			//[MWe]
		cr_out_solver.m_time_required_su = std::fmax(cr_out_solver.m_time_required_su, mc_gen3gas_receivers.at(i).ms_outputs.m_time_required_su);	//[s]
		cr_out_solver.m_W_dot_col_tracking += mc_pt_heliostatfields.at(i).ms_outputs.m_pparasi;

		q_dot_inc += mc_gen3gas_receivers.at(i).ms_outputs.m_q_dot_rec_inc;
		eta_thermal += mc_gen3gas_receivers.at(i).ms_outputs.m_eta_therm*mc_gen3gas_receivers.at(i).ms_outputs.m_Q_thermal;
		q_dot_thermal += mc_gen3gas_receivers.at(i).ms_outputs.m_Q_thermal;
		q_dot_startup += mc_gen3gas_receivers.at(i).ms_outputs.m_q_startup / (mc_gen3gas_receivers.at(i).ms_outputs.m_time_required_su / 3600.0);
		T_htf_out += mc_gen3gas_receivers.at(i).ms_outputs.m_T_salt_hot*norm;
		q_dot_pipe_loss += mc_gen3gas_receivers.at(i).ms_outputs.m_q_dot_piping_loss;
		q_dot_loss += mc_gen3gas_receivers.at(i).ms_outputs.m_q_rad_sum + mc_gen3gas_receivers.at(i).ms_outputs.m_q_conv_sum;
		q_dot_field_inc += mc_pt_heliostatfields.at(i).ms_outputs.m_q_dot_field_inc;
		eta_field += mc_pt_heliostatfields.at(i).ms_outputs.m_eta_field * mc_pt_heliostatfields.at(i).ms_outputs.m_q_dot_field_inc;
		sf_adjust_out += mc_pt_heliostatfields.at(i).ms_outputs.m_sf_adjust_out * mc_pt_heliostatfields.at(i).ms_outputs.m_q_dot_field_inc;
	}

	//The "total" mass flow is the max flow through the system, and this is governed by the first (north) receiver
	m_dot_htf = cr_out_solver.m_m_dot_salt_tot = mc_gen3gas_receivers.front().ms_outputs.m_m_dot_tot;		//[kg/hr]

	mc_reported_outputs.value(E_FIELD_Q_DOT_INC, q_dot_field_inc);	//[MWt]
	mc_reported_outputs.value(E_FIELD_ETA_OPT, eta_field);			//[-]
	mc_reported_outputs.value(E_FIELD_ADJUST, sf_adjust_out);			//[-]

	mc_reported_outputs.value(E_Q_DOT_INC, q_dot_inc);	//[MWt]
	mc_reported_outputs.value(E_ETA_THERMAL, eta_thermal);		//[-]
	mc_reported_outputs.value(E_Q_DOT_THERMAL, q_dot_thermal);	//[MWt]
	mc_reported_outputs.value(E_M_DOT_HTF, m_dot_htf);	//[kg/hr]
		// If startup, then timestep may have changed (why not report this from 222 in MWt?)
	mc_reported_outputs.value(E_Q_DOT_STARTUP, q_dot_startup);		//[MWt])
	mc_reported_outputs.value(E_T_HTF_IN, htf_state_in.m_temp);									//[C]
	mc_reported_outputs.value(E_T_HTF_OUT, T_htf_out);		//[C]
	mc_reported_outputs.value(E_Q_DOT_PIPE_LOSS, q_dot_pipe_loss);	//[MWt]
	mc_reported_outputs.value(E_Q_DOT_LOSS, q_dot_loss); //MWt
	return;
}

void C_csp_gas_collector_receiver::startup(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	C_csp_collector_receiver::S_csp_cr_out_solver &cr_out_solver,
	const C_csp_solver_sim_info &sim_info)
{
	// For now, define startup(...) shell that calls call() with operation mode defined.
	// Should eventually develop a startup method for the collector receiver

	// Set heliostat field call() parameters and solve
	C_csp_collector_receiver::S_csp_cr_inputs inputs;
	inputs.m_input_operation_mode = C_csp_collector_receiver::STARTUP;
	inputs.m_field_control = 1.0;

	call(weather, htf_state_in, inputs, cr_out_solver, sim_info);
}

void C_csp_gas_collector_receiver::on(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	double field_control,
	C_csp_collector_receiver::S_csp_cr_out_solver &cr_out_solver,
	//C_csp_collector_receiver::S_csp_cr_out_report &cr_out_report,
	const C_csp_solver_sim_info &sim_info)
{
	// For now, define on(...) shell that calls call() with operation mode defined.
	// Should eventually develop an 'on' method for the MSPT

	// Define 'C_csp_cr_inputs' for call(...)
	C_csp_collector_receiver::S_csp_cr_inputs inputs;
	inputs.m_input_operation_mode = C_csp_collector_receiver::ON;
	inputs.m_field_control = field_control;

	call(weather, htf_state_in, inputs, cr_out_solver, sim_info);
}

void C_csp_gas_collector_receiver::estimates(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	C_csp_collector_receiver::S_csp_cr_est_out &est_out,
	const C_csp_solver_sim_info &sim_info)
{
	// For now, define estimates(...) shell that calls call() with operation mode defined.
	// Should eventually develop an estimate(...) method for the MSPT
	
	C_csp_collector_receiver::S_csp_cr_inputs inputs;
	inputs.m_input_operation_mode = C_csp_collector_receiver::STEADY_STATE;
	inputs.m_field_control = 1.0;

	C_csp_collector_receiver::S_csp_cr_out_solver cr_out_solver;

	call(weather, htf_state_in, inputs, cr_out_solver, sim_info);

	int mode = get_operating_state();

	if( mode == C_csp_collector_receiver::ON )
	{
		est_out.m_q_dot_avail = cr_out_solver.m_q_thermal;			//[MWt]
		est_out.m_m_dot_avail = cr_out_solver.m_m_dot_salt_tot;		//[kg/hr]
		est_out.m_T_htf_hot = cr_out_solver.m_T_salt_hot;			//[C]
		est_out.m_q_startup_avail = 0.0;
	}
	else
	{
		est_out.m_q_startup_avail = cr_out_solver.m_q_thermal;		//[MWt]
		est_out.m_q_dot_avail = 0.0;
		est_out.m_m_dot_avail = 0.0;
		est_out.m_T_htf_hot = 0.0;
	}
}

double C_csp_gas_collector_receiver::calculate_optical_efficiency( const C_csp_weatherreader::S_outputs &weather, const C_csp_solver_sim_info &sim )
{
    /*
    Evaluate optical efficiency. This is a required function for the parent class, 
    but doesn't do much other than simply call the optical efficiency model in this case.
    */
	
	double eta_field_ave = 0.;

	for (size_t i = 0; i < mc_pt_heliostatfields.size(); i++)
	{
		mc_pt_heliostatfields.at(i).call(weather, 1., sim );
		eta_field_ave += mc_pt_heliostatfields.at(i).ms_outputs.m_eta_field;
	}

	return eta_field_ave / (double)mc_pt_heliostatfields.size();
}

double C_csp_gas_collector_receiver::get_collector_area()
{
	double A_sf_tot = 0.;
	for (size_t i = 0; i < mc_pt_heliostatfields.size(); i++)
		A_sf_tot += mc_pt_heliostatfields.at(i).ms_params.m_A_sf;

	return A_sf_tot;
}

double C_csp_gas_collector_receiver::calculate_thermal_efficiency_approx( const C_csp_weatherreader::S_outputs &weather, double q_inc )
{
    /* 
    A very approximate thermal efficiency used for quick optimization performance projections
    */

	return 0.85;

    //double T_eff = (mc_gen3gas_receivers.front().m_T_cold_des + mc_gen3gas_receivers.front().m_T_hot_des)*.55;

    //double T_amb = weather.m_tdry + 273.15;
    //double T_eff4 = T_eff * T_eff;
    //T_eff4 *= T_eff4;
    //double T_amb4 = T_amb * T_amb;
    //T_amb4 *= T_amb4;

    //double Arec = mc_gen3gas_receiver.m_d_rec * 3.1415 * mc_gen3gas_receiver.m_h_rec;

    //double q_rad = 5.67e-8*mc_gen3gas_receiver.m_epsilon * Arec * (T_eff4 - T_amb4) * 1.e-6;   //MWt

    //double v = weather.m_wspd;
    //double v2 = v*v;
    //double v3 = v2*v;

    //double q_conv = q_rad/2. * (-0.001129*v3 + 0.031229*v2 - 0.01822*v +0.962476);  //convection is about half radiation, scale by wind speed. surrogate regression from molten salt run.

    //return max(1. - (q_rad + q_conv)/q_inc, 0.);

}


void C_csp_gas_collector_receiver::converged()
{
	for (size_t i = 0; i < mc_gen3gas_receivers.size(); i++)
	{
		mc_pt_heliostatfields.at(i).converged();
		mc_gen3gas_receivers.at(i).converged();
	}

	// Hardcode to test...
	//mc_reported_outputs.set_timestep_output(E_Q_DOT_THERMAL, mc_gen3gas_receiver.ms_outputs.m_Q_thermal);	//[MWt]
	mc_reported_outputs.set_timestep_outputs();
}

void C_csp_gas_collector_receiver::write_output_intervals(double report_time_start,
	const std::vector<double> & v_temp_ts_time_end, double report_time_end)
{
	mc_reported_outputs.send_to_reporting_ts_array(report_time_start,
		v_temp_ts_time_end, report_time_end);
}