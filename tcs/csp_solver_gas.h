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

#ifndef __csp_solver_gas_
#define __csp_solver_gas_

#include <numeric>
#include <limits>
#include <memory>

#include "lib_weatherfile.h"
#include "csp_solver_util.h"
#include "csp_solver_core.h"

#include "numeric_solvers.h"


class C_csp_solver_gas
{

public:
	
	class C_solver_outputs
	{
	public:
		enum
		{
			// Ouputs that are NOT reported as weighted averages
				// Simulation
			TIME_FINAL,       //[hr] Simulation timestep
				// Weather Reader
			MONTH,            //[-] Month of year
			HOUR_DAY,         //[hr] hour of day
				// Controller
			ERR_M_DOT,        //[-] Relative mass conservation error
			ERR_Q_DOT,        //[-] Relative energy conservation error
			N_OP_MODES,       //[-] Number of subtimesteps in reporting timestep
			OP_MODE_1,        //[-] First operating mode in reporting timestep - always should be valid
			OP_MODE_2,        //[-] 2nd operating mode in reporting timestep - not always populated
			OP_MODE_3,        //[-] 3rd operating mode in reporting timestep - usually NOT populated

			// **************************************************************
			//      ONLY instantaneous outputs that are reported as the first value
			//        if multiple csp-timesteps for one reporting timestep
			// **************************************************************
			TOU_PERIOD,                 //[-] CSP operating TOU period
			PRICING_MULT,               //[-] PPA price multiplier
			PC_Q_DOT_SB,                //[MWt] PC required standby thermal power
			PC_Q_DOT_MIN,               //[MWt] PC required min thermal power
			PC_Q_DOT_TARGET,            //[MWt] PC target thermal power
			PC_Q_DOT_MAX,               //[MWt] PC allowable max thermal power
			CTRL_IS_REC_SU,             //[-] Control decision: is receiver startup allowed?
			CTRL_IS_PC_SU,              //[-] Control decision: is power cycle startup allowed?
			CTRL_IS_PC_SB,              //[-] Control decision: is power cycle standby allowed?
			EST_Q_DOT_CR_SU,            //[MWt] Estimate receiver startup thermal power
			EST_Q_DOT_CR_ON,            //[MWt] Estimate receiver thermal power to HTF
			EST_Q_DOT_DC,               //[MWt] Estimate max TES dc thermal power
			EST_Q_DOT_CH,               //[MWt] Estimate max TES ch thermal power
			CTRL_OP_MODE_SEQ_A,         //[-] First 3 operating modes tried
			CTRL_OP_MODE_SEQ_B,         //[-] Next 3 operating modes tried
			CTRL_OP_MODE_SEQ_C,         //[-] Final 3 operating modes tried

			// **************************************************************
			//      Outputs that are reported as weighted averages if 
			//       multiple csp-timesteps for one reporting timestep
			// **************************************************************
				// Weather Reader			
			SOLZEN,           //[deg] Solar zenith angle
			SOLAZ,            //[deg] Solar azimuth angle
			BEAM,			  //[W/m^2] Resource beam normal irradiance
			TDRY,             //[C] Dry bulb temperature
			TWET,             //[C] Wet bulb temperature
			RH,	              //[-] Relative humidity
			WSPD,             //[m/s] Wind speed
			PRES,             //[mbar] Atmospheric pressure
				// Controller and Storage
			CR_DEFOCUS,       //[-] Field optical focus fraction
			TES_Q_DOT_LOSS,       //[MWt] TES thermal losses
			TES_W_DOT_HEATER,     //[MWe] TES freeze protection power
			TES_T_HOT,            //[C] TES final hot tank temperature
			TES_T_COLD,           //[C] TES final cold tank temperature
			TES_Q_DOT_DC,         //[MWt] TES discharge thermal power
			TES_Q_DOT_CH,         //[MWt] TES charge thermal power
			TES_E_CH_STATE,       //[MWht] TES charge state at the end of the time step
			TES_M_DOT_DC,         //[MWt] TES discharge mass flow rate
			TES_M_DOT_CH,         //[MWt] TES charge mass flow rate
			COL_W_DOT_TRACK,      //[MWe] Parasitic collector tracking, startup, stow power consumption
			CR_W_DOT_PUMP,        //[MWe] Parasitic tower HTF pump power
			SYS_W_DOT_PUMP,       //[MWe] Parasitic PC and TES HTF pump power
			PC_W_DOT_COOLING,     //[MWe] Parasitic condenser operation power
			SYS_W_DOT_FIXED,      //[MWe] Parasitic fixed power consumption
			SYS_W_DOT_BOP,        //[MWe] Parasitic BOP power consumption
			W_DOT_NET             //[MWe] System total electric power to grid
		};
	};
	
	C_csp_reported_outputs mc_reported_outputs;

	struct S_sim_setup
	{
		double m_sim_time_start;	//[s]
		double m_sim_time_end;		//[s]
		double m_report_step;		//[s]

		S_sim_setup()
		{
			m_sim_time_start = m_sim_time_end = m_report_step = std::numeric_limits<double>::quiet_NaN();
		}
	};

	class C_csp_solver_kernel
	{
	private:
		S_sim_setup ms_sim_setup;

		C_timestep_fixed mc_ts_weatherfile;

		C_timestep_fixed mc_ts_sim_baseline;

	public:
			
		C_csp_solver_sim_info mc_sim_info;

		void init(C_csp_solver_gas::S_sim_setup & sim_setup, double wf_step /*s*/, 
			double baseline_step /*s*/, C_csp_messages & csp_messages);

		void wf_step_forward();

		void baseline_step_forward();

		double get_wf_end_time();
		double get_wf_step();

		double get_baseline_end_time();
		double get_baseline_step();		

		const S_sim_setup * get_sim_setup();

	};
	
	struct S_csp_system_params
	{
		double m_pb_fixed_par;		//[MWe/MWcap]
		
		double m_bop_par;			//[MWe/MWcap]
		double m_bop_par_f;			//[-]
		double m_bop_par_0;			//[-]
		double m_bop_par_1;			//[-]
		double m_bop_par_2;			//[-]

		S_csp_system_params()
		{
			m_pb_fixed_par =

			m_bop_par = m_bop_par_f = m_bop_par_0 = m_bop_par_1 = m_bop_par_2 = std::numeric_limits<double>::quiet_NaN();
		}
	};

private:
	C_csp_weatherreader &mc_weather;
	C_csp_collector_receiver &mc_collector_receiver;
	C_csp_power_cycle &mc_power_cycle;
	C_csp_tes &mc_tes;
	C_csp_tou &mc_tou;

	S_csp_system_params & ms_system_params;

	C_csp_solver_htf_1state mc_cr_htf_state_in;
	C_csp_collector_receiver::S_csp_cr_out_solver mc_cr_out_solver;

	C_csp_solver_htf_1state mc_pc_htf_state_in;
	C_csp_power_cycle::S_control_inputs mc_pc_inputs;
	C_csp_power_cycle::S_csp_pc_out_solver mc_pc_out_solver;

	C_csp_solver_htf_state mc_tes_ch_htf_state;
	C_csp_solver_htf_state mc_tes_dc_htf_state;
	C_csp_tes::S_csp_tes_outputs mc_tes_outputs;

    C_csp_tou::S_csp_tou_outputs mc_tou_outputs;

	C_csp_solver_gas::C_csp_solver_kernel mc_kernel;

	// Hierarchy logic
	bool m_is_CR_SU__PC_OFF__TES_OFF_avail;
	bool m_is_CR_SU__PC_MIN__TES_EMPTY_avail;
	bool m_is_CR_SU__PC_SB__TES_DC_avail;
	bool m_is_CR_SU__PC_TARGET__TES_DC_avail;
	bool m_is_CR_SU__PC_RM_LO__TES_EMPTY_avail;
	bool m_is_CR_SU__PC_SU__TES_DC_avail;
	
	bool m_is_CR_ON__PC_OFF__TES_CH_avail;
	bool m_is_CR_ON__PC_SU__TES_OFF_avail;
	bool m_is_CR_ON__PC_SU__TES_CH_avail;
	bool m_is_CR_ON__PC_MIN__TES_EMPTY_avail;
	bool m_is_CR_ON__PC_SB__TES_DC_avail;
	bool m_is_CR_ON__PC_SB__TES_OFF_avail;
	bool m_is_CR_ON__PC_SB__TES_CH_avail;
	bool m_is_CR_ON__PC_SB__TES_FULL_avail;
	bool m_is_CR_ON__PC_RM_HI__TES_FULL_avail;
	bool m_is_CR_ON__PC_RM_HI__TES_OFF_avail_HI_SIDE;
	bool m_is_CR_ON__PC_RM_HI__TES_OFF_avail_LO_SIDE;
	bool m_is_CR_ON__PC_RM_LO__TES_OFF_avail;
	bool m_is_CR_ON__PC_RM_LO__TES_EMPTY_avail;
	bool m_is_CR_ON__PC_TARGET__TES_CH_avail_HI_SIDE;
	bool m_is_CR_ON__PC_TARGET__TES_CH_avail_LO_SIDE;
	bool m_is_CR_ON__PC_TARGET__TES_DC_avail;

	bool m_is_CR_DF__PC_OFF__TES_FULL_avail;
	bool m_is_CR_DF__PC_SU__TES_FULL_avail;
	bool m_is_CR_DF__PC_SU__TES_OFF_avail;
	bool m_is_CR_DF__PC_MAX__TES_OFF_avail;
	bool m_is_CR_DF__PC_MAX__TES_FULL_avail;
	
	bool m_is_CR_OFF__PC_SU__TES_DC_avail;
	bool m_is_CR_OFF__PC_MIN__TES_EMPTY_avail;
	bool m_is_CR_OFF__PC_RM_LO__TES_EMPTY_avail;
	bool m_is_CR_OFF__PC_TARGET__TES_DC_avail;
	bool m_is_CR_OFF__PC_SB__TES_DC_avail;


	// member string for exception messages
	std::string error_msg;

		// Collector receiver design parameters
	double m_T_htf_cold_des;			//[K]
	double m_P_cold_des;				//[kPa]
	double m_x_cold_des;				//[-]
	double m_q_dot_rec_des;				//[MW]
	double m_A_aperture;				//[m2]

		// Power cycle design parameters
	double m_cycle_W_dot_des;			//[MW]
	double m_cycle_eta_des;				//[-]
	double m_cycle_q_dot_des;			//[MW]
	double m_cycle_max_frac;			//[-]
	double m_cycle_cutoff_frac;			//[-]
	double m_cycle_sb_frac_des;			//[-]
	double m_cycle_T_htf_hot_des;		//[K]
	double m_cycle_P_hot_des;			//[kPa]
	double m_cycle_x_hot_des;			//[-]
	double m_m_dot_pc_des;				//[kg/hr]
	double m_m_dot_pc_min;				//[kg/hr]
	double m_m_dot_pc_max;				//[kg/hr]

		// Storage logic
	bool m_is_tes;			//[-] True: plant has storage

		// Reporting and Output Tracking
	int m_i_reporting;					//[-]
	double m_report_time_start;			//[s]
	double m_report_time_end;			//[s]
	double m_report_step;				//[s]

		// Estimates to use
	double m_T_htf_pc_cold_est;			//[C]

		// Reset hierarchy logic
	void reset_hierarchy_logic();
	void turn_off_plant();

	void solver_cr_to_pc_to_cr(int pc_mode, double field_control_in, double tol, int &exit_mode, double &exit_tolerance);
	 
	void solver_pc_su_controlled__tes_dc(double step_tol /*s*/,
		double &time_pc_su /*s*/, 
		int & exit_mode, double &T_pc_in_exit_tolerance);

	int solver_pc_fixed__tes_empty(double q_dot_pc_fixed /*MWt*/,
		double tol,
		double & time_tes_dc);

	int solver_cr_on__pc_match__tes_full(int pc_mode, double defocus_in);

	enum E_solver_outcomes
	{
		CSP_NO_SOLUTION,	// Models did not provide enough information with which to iterate on T_rec_in
		POOR_CONVERGENCE,	// Models solved, but convergence on T_rec_in was not within specified tolerance
		CSP_CONVERGED,		// Models solved; convergence within specified tolerance
		KNOW_NEXT_MODE,		// Models did not solve, but failure mode allowed next mode to be determined
		UNDER_TARGET_PC,	// Models solved, but could not converge because the operating mode did not allow enough thermal power to go to power cycle
		OVER_TARGET_PC,		// Models solved, but could not converge because the operating mode could not reduce the mass flow rate enough to the power cycle
		REC_IS_OFF			// Collector-receiver model did not produce power
	};

	// Solved Controller Variables
	double m_defocus;		//[-] (1..0) Should only be less than 1 if receiver is on, but defocused
	
	std::vector<double> mv_time_local;

	bool(*mpf_callback)(std::string &log_msg, std::string &progress_msg, void *data, double progress, int log_type);
	void *mp_cmod_active;

	void send_callback(double percent);

public:

	// Class to save messages for up stream classes
	C_csp_messages mc_csp_messages;

	// Vector to track operating modes
	std::vector<int> m_op_mode_tracking;

	enum tech_operating_modes
	{
		ENTRY_MODE = 0,
		
		CR_OFF__PC_OFF__TES_OFF,
		CR_SU__PC_OFF__TES_OFF,
		CR_ON__PC_SU__TES_OFF,
		CR_ON__PC_SB__TES_OFF,
		
		CR_ON__PC_RM_HI__TES_OFF,
		CR_ON__PC_RM_LO__TES_OFF,
		
		CR_DF__PC_MAX__TES_OFF,

		CR_OFF__PC_SU__TES_DC,
		CR_ON__PC_OFF__TES_CH,

		SKIP_10,

		CR_ON__PC_TARGET__TES_CH,
		CR_ON__PC_TARGET__TES_DC,

		CR_ON__PC_RM_LO__TES_EMPTY,

		CR_DF__PC_OFF__TES_FULL,
		
		CR_OFF__PC_SB__TES_DC,
		CR_OFF__PC_MIN__TES_EMPTY,
		CR_OFF__PC_RM_LO__TES_EMPTY,

		CR_ON__PC_SB__TES_CH,
		CR_SU__PC_MIN__TES_EMPTY,

		SKIP_20,

		CR_SU__PC_SB__TES_DC,
		CR_ON__PC_SB__TES_DC,
		CR_OFF__PC_TARGET__TES_DC,
		CR_SU__PC_TARGET__TES_DC,
		CR_ON__PC_RM_HI__TES_FULL,

		CR_ON__PC_MIN__TES_EMPTY,

		CR_SU__PC_RM_LO__TES_EMPTY,

		CR_DF__PC_MAX__TES_FULL,

		CR_ON__PC_SB__TES_FULL,

		SKIP_30,

		CR_SU__PC_SU__TES_DC,

		CR_ON__PC_SU__TES_CH,

		CR_DF__PC_SU__TES_FULL,

		CR_DF__PC_SU__TES_OFF
	};

    static std::string tech_operating_modes_str[];
    
	C_csp_solver_gas(C_csp_weatherreader &weather,
		C_csp_collector_receiver &collector_receiver,
		C_csp_power_cycle &power_cycle,
		C_csp_tes &tes,
		C_csp_tou &tou,
		S_csp_system_params &system,
		bool(*pf_callback)(std::string &log_msg, std::string &progress_msg, void *data, double progress, int out_type) = 0,
		void *p_cmod_active = 0);

	~C_csp_solver_gas(){};

	void init();

	void Ssimulate(C_csp_solver_gas::S_sim_setup & sim_setup);

	int steps_per_hour();

	double get_cr_aperture_area();

	// Output vectors
	// Need to be sure these are always up-to-date as multiple operating modes are tested during one timestep
	std::vector< std::vector< double > > mvv_outputs_temp;

	// *****************************
	// *****************************
	// Solvers
	class C_MEQ_cr_on__pc_q_dot_max__tes_off__defocus : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		int m_pc_mode;				//[-]
		double m_q_dot_max;			//[MWt]

	public:
		C_MEQ_cr_on__pc_q_dot_max__tes_off__defocus(C_csp_solver_gas *pc_csp_solver, int pc_mode /*-*/, double q_dot_max /*MWt*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_pc_mode = pc_mode;
			m_q_dot_max = q_dot_max;			
		}

		virtual int operator()(double defocus /*-*/, double *q_dot_pc /*MWt*/);
	};

	class C_mono_eq_cr_to_pc_to_cr : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		int m_pc_mode;					//[-]
		double m_P_field_in;			//[kPa]
		double m_x_field_in;			//[-]
		double m_field_control_in;		//[-]

	public:
		C_mono_eq_cr_to_pc_to_cr(C_csp_solver_gas *pc_csp_solver, int pc_mode /*-*/,
			double P_field_in /*kPa*/, double x_field_in /*-*/, double field_control_in /*-*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_pc_mode = pc_mode;
			m_P_field_in = P_field_in;
			m_x_field_in = x_field_in;
			m_field_control_in = field_control_in;
		}
		
		virtual int operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/);
	};

	class C_mono_eq_pc_su_cont_tes_dc : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;

	public:
		C_mono_eq_pc_su_cont_tes_dc(C_csp_solver_gas *pc_csp_solver)
		{
			mpc_csp_solver = pc_csp_solver;
			m_time_pc_su = std::numeric_limits<double>::quiet_NaN();
		}

		double m_time_pc_su;		//[s] power cycle model returns MIN(time required to completely startup, full timestep duration)

		virtual int operator()(double T_htf_hot /*C*/, double *diff_T_htf_hot /*-*/);
	};

	class C_mono_eq_pc_target_tes_dc__m_dot : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		double m_T_htf_cold;		//[C]
		int m_pc_mode;				//[-]

	public:
		C_mono_eq_pc_target_tes_dc__m_dot(C_csp_solver_gas *pc_csp_solver, 
				int pc_mode /*-*/, double T_htf_cold /*C*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_T_htf_cold = T_htf_cold;
			m_pc_mode = pc_mode;
		}

		virtual int operator()(double m_dot_htf /*kg/hr*/, double *q_dot_pc /*MWt*/);
	};

	class C_mono_eq_pc_target_tes_dc__T_cold : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		int m_pc_mode;			//[-]
		double m_q_dot_target;		//[MWt]

	public:
		C_mono_eq_pc_target_tes_dc__T_cold(C_csp_solver_gas *pc_csp_solver, 
				int pc_mode /*-*/, double q_dot_target /*MWt*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_pc_mode = pc_mode;	//[-]
			m_q_dot_target = q_dot_target;		//[MWt]
			m_q_dot_calc = std::numeric_limits<double>::quiet_NaN();	//[MWt]
			m_m_dot_calc = std::numeric_limits<double>::quiet_NaN();	//[kg/hr]
		}

		double m_q_dot_calc;		//[MWt]
		double m_m_dot_calc;		//[kg/hr]

		virtual int operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/);
	};

	class C_mono_eq_pc_match_tes_empty : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;

	public:
		C_mono_eq_pc_match_tes_empty(C_csp_solver_gas *pc_csp_solver)
		{
			mpc_csp_solver = pc_csp_solver;
		}

		virtual int operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/);
	};

	class C_mono_eq_cr_on_pc_su_tes_ch : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;

	public:
		C_mono_eq_cr_on_pc_su_tes_ch(C_csp_solver_gas *pc_csp_solver)
		{
			mpc_csp_solver = pc_csp_solver;
			m_step_pc_su = std::numeric_limits<double>::quiet_NaN();
		}

		double m_step_pc_su;	//[s]

		virtual int operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/);
	};

	class C_mono_eq_pc_target__m_dot : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		int m_pc_mode;		//[-]
		double m_T_htf_hot;	//[C]

	public:
		C_mono_eq_pc_target__m_dot(C_csp_solver_gas *pc_csp_solver, int pc_mode, double T_htf_hot /*C*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_pc_mode = pc_mode;
			m_T_htf_hot = T_htf_hot;	//[C]
		}

		virtual int operator()(double m_dot_htf_pc /*kg/hr*/, double *q_dot_pc /*MWt*/);
	};

	class C_mono_eq_cr_on_pc_target_tes_ch__T_cold : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		int m_pc_mode;			//[-]
		double m_q_dot_target;		//[MWt]
		double m_defocus;			//[-]

	public:
		C_mono_eq_cr_on_pc_target_tes_ch__T_cold(C_csp_solver_gas *pc_csp_solver, 
									int pc_mode, double q_dot_target /*MWt*/, double defocus /*-*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_pc_mode = pc_mode;				//[-]
			m_q_dot_target = q_dot_target;		//[MWt]
			m_defocus = defocus;				//[-]
		}

		virtual int operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/);
	};

	class C_mono_eq_cr_on_pc_match_tes_empty : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		double m_defocus;		//[-]
		
	public:
		C_mono_eq_cr_on_pc_match_tes_empty(C_csp_solver_gas *pc_csp_solver, double defocus /*-*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_defocus = defocus;	//[-]
		}

		virtual int operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/);
	};

	class C_mono_eq_cr_on_pc_target_tes_dc : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		int m_pc_mode;			//[-]
		double m_defocus;		//[-]
		double m_q_dot_target;	//[MWt]

	public:
		C_mono_eq_cr_on_pc_target_tes_dc(C_csp_solver_gas *pc_csp_solver,
			int pc_mode, double q_dot_target /*MWt*/, double defocus /*-*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_pc_mode = pc_mode;				//[-]
			m_defocus = defocus;				//[-]
			m_q_dot_target = q_dot_target;		//[MWt]
		}

		virtual int operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/);
	};

	class C_mono_eq_pc_target__m_dot_fixed_plus_tes_dc : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		int m_pc_mode;			//[-]
		double m_T_htf_cold;	//[C]
		double m_T_htf_fixed_hot;	//[C]
		double m_m_dot_htf_fixed;	//[kg/hr]

	public:
		C_mono_eq_pc_target__m_dot_fixed_plus_tes_dc(C_csp_solver_gas *pc_csp_solver,
						int pc_mode, double T_htf_cold /*C*/,
						double T_htf_fixed_hot /*C*/, double m_dot_htf_fixed /*kg/hr*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_pc_mode = pc_mode;
			m_T_htf_cold = T_htf_cold;				//[C]
			m_T_htf_fixed_hot = T_htf_fixed_hot;	//[C]
			m_m_dot_htf_fixed = m_dot_htf_fixed;	//[kg/hr]
		}

		virtual int operator()(double m_dot_tes_dc /*kg/hr*/, double *q_dot_pc);
	};

	class C_mono_eq_pc_target_tes_empty__x_step : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		double m_T_htf_cold;			//[C]

	public:
		C_mono_eq_pc_target_tes_empty__x_step(C_csp_solver_gas *pc_csp_solver,
									double T_htf_cold /*C*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_T_htf_cold = T_htf_cold;				//[C]
		}

		virtual int operator()(double step /*s*/, double *q_dot_pc /*MWt*/);
	};

	class C_mono_eq_pc_target_tes_empty__T_cold : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		double m_q_dot_pc_target;		//[MWt]

	public:
		C_mono_eq_pc_target_tes_empty__T_cold(C_csp_solver_gas *pc_csp_solver,
								double q_dot_pc_target /*MWt*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_q_dot_pc_target = q_dot_pc_target;	//[MWt]
			m_step = std::numeric_limits<double>::quiet_NaN();
		}

		double m_step;		//[s]

		virtual int operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/);

		void solve_pc(double step /*s*/);
	};

	class C_mono_eq_cr_on__pc_match__tes_full : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		int m_pc_mode;					//[-]
		double m_defocus;				//[-]

	public:
		C_mono_eq_cr_on__pc_match__tes_full(C_csp_solver_gas *pc_csp_solver,
										int pc_mode, double defocus /*-*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_pc_mode = pc_mode;
			m_defocus = defocus;
		}

		virtual int operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/);
	};

	class C_mono_eq_cr_on__pc_max_m_dot__tes_full : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		int m_pc_mode;					//[-]
		double m_defocus;				//[-]

	public:
		C_mono_eq_cr_on__pc_max_m_dot__tes_full(C_csp_solver_gas *pc_csp_solver,
									int pc_mode, double defocus /*-*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_pc_mode = pc_mode;
			m_defocus = defocus;
		}

		virtual int operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/);
	};

	class C_mono_eq_cr_on__pc_target__tes_full__defocus : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		int m_pc_mode;			//[-]
		double m_q_dot_max;		//[MWt]

	public:
		C_mono_eq_cr_on__pc_target__tes_full__defocus(C_csp_solver_gas *pc_csp_solver,
			int pc_mode, double q_dot_max /*MWt*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_pc_mode = pc_mode;
			m_q_dot_max = q_dot_max;		//[MWt]
		}

		virtual int operator()(double defocus /*-*/, double *q_dot_pc /*MWt*/);
	};

	class C_mono_eq_cr_on__pc_m_dot_max__tes_full_defocus : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		int m_pc_mode;		//[-]

	public:
		C_mono_eq_cr_on__pc_m_dot_max__tes_full_defocus(C_csp_solver_gas *pc_csp_solver,
			int pc_mode)
		{
			mpc_csp_solver = pc_csp_solver;
			m_pc_mode = pc_mode;
		}

		virtual int operator()(double defocus /*-*/, double *m_dot_bal /*-*/);
	};

	class C_mono_eq_cr_on__pc_match_m_dot_ceil__tes_full : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		int m_pc_mode;					//[-]
		double m_defocus;				//[-]

	public:
		C_mono_eq_cr_on__pc_match_m_dot_ceil__tes_full(C_csp_solver_gas *pc_csp_solver,
			int pc_mode, double defocus /*-*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_pc_mode = pc_mode;
			m_defocus = defocus;
		}

		virtual int operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/);
	};

	class C_MEQ_cr_on__pc_m_dot_max__tes_off__defocus : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		int m_pc_mode;		//[-]

	public:
		C_MEQ_cr_on__pc_m_dot_max__tes_off__defocus(C_csp_solver_gas *pc_csp_solver,
			int pc_mode)
		{
			mpc_csp_solver = pc_csp_solver;
			m_pc_mode = pc_mode;
		}

		virtual int operator()(double defocus /*-*/, double *m_dot_bal /*-*/);
	};

	class C_MEQ_cr_on__pc_max_m_dot__tes_off__T_htf_cold : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		int m_pc_mode;					//[-]
		double m_defocus;				//[-]

	public:
		C_MEQ_cr_on__pc_max_m_dot__tes_off__T_htf_cold(C_csp_solver_gas *pc_csp_solver,
			int pc_mode, double defocus /*-*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_pc_mode = pc_mode;
			m_defocus = defocus;
		}

		virtual int operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/);
	};

	class C_MEQ_cr_on__pc_off__tes_ch__T_htf_cold : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		double m_defocus;

	public:
		C_MEQ_cr_on__pc_off__tes_ch__T_htf_cold(C_csp_solver_gas *pc_csp_solver,
							double defocus /*-*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_defocus = defocus;				//[-]
		}

		virtual int operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/);
	};

	class C_MEQ_cr_on__pc_target__tes_empty__T_htf_cold : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		double m_defocus;
		double m_q_dot_pc_target;	//[MWt]

	public:
		C_MEQ_cr_on__pc_target__tes_empty__T_htf_cold(C_csp_solver_gas *pc_csp_solver,
			double defocus /*-*/, double q_dot_pc_target /*MWt*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_defocus = defocus;
			m_q_dot_pc_target = q_dot_pc_target;
			m_step = std::numeric_limits<double>::quiet_NaN();
		}

		double m_step;	//[s]

		virtual int operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/);

		void solve_pc(double step /*s*/, double T_htf_pc_hot /*C*/, double m_dot_htf_pc /*kg/hr*/);
	};

	class C_MEQ_cr_on__pc_target__tes_empty__step : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		double m_defocus;				//[-]
		double m_T_htf_cold;			//[C]

	public:
		C_MEQ_cr_on__pc_target__tes_empty__step(C_csp_solver_gas *pc_csp_solver,
			double defocus /*-*/, double T_htf_cold /*C*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_defocus = defocus;
			m_T_htf_cold = T_htf_cold;				//[C]
			m_m_dot_pc = std::numeric_limits<double>::quiet_NaN();		//[kg/hr]
			m_T_htf_pc_hot = std::numeric_limits<double>::quiet_NaN();	//[MWt]
		}

		double m_m_dot_pc;			//[kg/hr]
		double m_T_htf_pc_hot;		//[C]

		virtual int operator()(double step /*s*/, double *q_dot_pc /*MWt*/);
	};

	class C_MEQ_cr_df__pc_off__tes_full__T_cold : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;
		double m_defocus;		//[-]

	public:
		C_MEQ_cr_df__pc_off__tes_full__T_cold(C_csp_solver_gas *pc_csp_solver,
			double defocus /*-*/)
		{
			mpc_csp_solver = pc_csp_solver;
			m_defocus = defocus;
		}

		virtual int operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/);
	};

	class C_MEQ_cr_df__pc_off__tes_full__defocus : public C_monotonic_equation
	{
	private:
		C_csp_solver_gas *mpc_csp_solver;

	public:
		C_MEQ_cr_df__pc_off__tes_full__defocus(C_csp_solver_gas *pc_csp_solver)
		{
			mpc_csp_solver = pc_csp_solver;
		}

		virtual int operator()(double defocus /*-*/, double *diff_m_dot /*-*/);
	};

};


#endif //
