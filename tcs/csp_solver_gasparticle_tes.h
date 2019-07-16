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

#ifndef __csp_solver_gasparticle_tes_
#define __csp_solver_gasparticle_tes_

#include "csp_solver_core.h"
#include "csp_solver_util.h"
#include "csp_solver_two_tank_tes.h"        //for C_heat_exchanger and C_storage_tank classes

#include "sam_csp_util.h"

class C_csp_gasparticle_tes : public C_csp_tes
{
private:

	HTFProperties mc_gas_htfProps;		// Instance of HTFProperties class for gas HTF
	HTFProperties mc_particle_htfProps;		// Instance of HTFProperties class for storage media
	
    std::vector<C_heat_exchanger> mc_charge_hxs;
    C_heat_exchanger mc_discharge_hx;

	C_storage_tank mc_cold_tank;			// Instance of storage tank class for the cold tank
	C_storage_tank mc_hot_tank;				// Instance of storage tank class for the hot tank	

	// member string for exception messages
	std::string error_msg;

	// Timestep data
	double m_m_dot_tes_dc_max;  //[kg/s] TES discharge available from the SYSTEM (field side of HX if there is one)
	double m_m_dot_tes_ch_max;  //[kg/s] TES charge that can be sent to the SYSTEM (field side of HX if there is one)

	// Member data
	double m_vol_tank;			//[m3] volume of *one temperature*, i.e. vol_tank = total cold storage = total hot storage
	double m_V_tank_active;		//[m^3] available volume (considering h_min) of *one temperature*
	double m_q_pb_design;		//[Wt] thermal power to power cycle at design
	double m_V_tank_hot_ini;	//[m^3] Initial volume in hot storage tank
    double m_d_silo;            //[m] Silo diameter (calculated)
    double m_h_silo;            //[m] Silo height (calculated)

    // Monotonic equation solvers
    class C_MEQ_indirect_tes_discharge : public C_monotonic_equation
    {
    private:
        C_csp_gasparticle_tes *mpc_csp_gasparticle_tes;
        double m_timestep;
        double m_T_amb;
        double m_T_cold_field;
        double m_m_dot_field;

    public:
        C_MEQ_indirect_tes_discharge(C_csp_gasparticle_tes *pc_csp_gasparticle_tes, double timestep, double T_amb,
            double T_cold_field, double m_dot_field)
        {
            mpc_csp_gasparticle_tes = pc_csp_gasparticle_tes;
            m_timestep = timestep;
            m_T_amb = T_amb;
            m_T_cold_field = T_cold_field;
            m_m_dot_field = m_dot_field;
        }

        virtual int operator()(double m_dot_tank /*kg/s*/, double *m_dot_bal /*-*/);
    };

    class C_MEQ_indirect_tes_charge : public C_monotonic_equation
    {
    private:
        C_csp_gasparticle_tes *mpc_csp_gasparticle_tes;
        double m_timestep;
        double m_T_amb;
        double m_T_hot_field;
        double m_m_dot_field;

    public:
        C_MEQ_indirect_tes_charge(C_csp_gasparticle_tes *pc_csp_gasparticle_tes, double timestep, double T_amb,
            double T_hot_field, double m_dot_field)
        {
            mpc_csp_gasparticle_tes = pc_csp_gasparticle_tes;
            m_timestep = timestep;
            m_T_amb = T_amb;
            m_T_hot_field = T_hot_field;
            m_m_dot_field = m_dot_field;
        }

        virtual int operator()(double m_dot_tank /*kg/s*/, double *m_dot_bal /*-*/);
    };

public:

	// Class to save messages for up stream classes
	C_csp_messages mc_csp_messages;

	struct S_params
	{
		int m_field_fl;
		util::matrix_t<double> m_field_fl_props;

		int m_tes_fl;
		util::matrix_t<double> m_tes_fl_props;

		double m_W_dot_pc_design;   //[MWe] Design point gross power cycle output
		double m_eta_pc;            //[-] Design point power cycle thermal efficiency
		double m_solarm;			//[-] solar multiple
		double m_ts_hours;			//[hr] hours of storage at design power cycle operation		
		double m_aspect_tank;			//[-] tank height / tank diameter
		double m_u_tank;			//[W/m^2-K]
        double m_dt_charge;			//[C] Approach temperature of charge heat exchangers
        double m_dt_discharge;			//[C] Approach temperature of discharge heat exchanger
		double m_T_cycle_out_des;	//[C] Cycle cold return temperature to discharge heat exchanger. convert to K in init()
		double m_T_receiver_out_des;	//[C] Receiver hot outlet temperature. convert to K in init()
		double m_T_tank_hot_ini;	//[C] Initial temperature in hot storage tank
		double m_T_tank_cold_ini;	//[C] Initial temperature in cold storage cold
		double m_h_tank_min;		//[m] Minimum allowable HTF height in storage tank
		double m_f_V_hot_ini;       //[%] Initial fraction of available volume that is hot
		double m_htf_pump_coef;		//[kW/kg/s] Pumping power to move 1 kg/s of HTF through power cycle
        double m_tes_pump_coef;		//[kW/kg/s] Pumping power to move 1 kg/s of HTF through tes loop
        double T_tank_hot_inlet_min; //[C] Minimum field htf temperature that may enter the hot tank
        int m_tank_pairs;			//[-]

		S_params()
		{
			m_field_fl = m_tes_fl = m_tank_pairs = -1;		

			m_W_dot_pc_design = m_eta_pc = m_solarm = m_ts_hours = m_aspect_tank = m_u_tank = m_dt_charge = m_dt_discharge = m_T_cycle_out_des = m_T_receiver_out_des = m_T_tank_hot_ini =
				m_T_tank_cold_ini = m_h_tank_min = m_f_V_hot_ini = m_htf_pump_coef = m_tes_pump_coef = T_tank_hot_inlet_min = std::numeric_limits<double>::quiet_NaN();
		}
	};

	S_params ms_params;

	C_csp_gasparticle_tes();

	~C_csp_gasparticle_tes(){};

	virtual void init();

	virtual bool does_tes_exist();

	virtual double get_hot_temp();

	virtual double get_cold_temp();



    virtual double get_initial_charge_energy(); //MWh

    virtual double get_min_charge_energy(); //MWh

    virtual double get_max_charge_energy(); //MWh

    virtual double get_degradation_rate();  // s^-1

	virtual void discharge_avail_est(double T_cold_K, double step_s, double &q_dot_dc_est, double &m_dot_field_est, double &T_hot_field_est);

	virtual void charge_avail_est(double T_hot_K, double step_s, double &q_dot_ch_est, double &m_dot_field_est, double &T_cold_field_est);

	// Calculate pumping power...???
	virtual bool discharge(double timestep /*s*/, double T_amb /*K*/, double m_dot_htf_in /*kg/s*/, double T_htf_cold_in, double & T_htf_hot_out /*K*/, C_csp_tes::S_csp_tes_outputs &outputs);

	virtual void discharge_full(double timestep /*s*/, double T_amb /*K*/, double T_htf_cold_in, double & T_htf_hot_out /*K*/, double & m_dot_htf_out /*kg/s*/, C_csp_tes::S_csp_tes_outputs &outputs);

	virtual bool charge(C_heat_exchanger& charge_hx, double timestep /*s*/, double T_amb /*K*/, double m_dot_htf_in /*kg/s*/, double T_htf_hot_in, double & T_htf_cold_out /*K*/, C_csp_tes::S_csp_tes_outputs &outputs);

	virtual void charge_full(double timestep /*s*/, double T_amb /*K*/, double T_htf_hot_in /*K*/, double & T_htf_cold_out /*K*/, double & m_dot_htf_out /*kg/s*/, C_csp_tes::S_csp_tes_outputs &outputs);
	
	virtual void idle(double timestep, double T_amb, C_csp_tes::S_csp_tes_outputs &outputs);
	
	virtual void converged();
};

void gasparticle_tes_sizing(HTFProperties &tes_htf_props, double Q_tes_des /*MWt-hr*/, double T_tes_hot /*K*/,
		double T_tes_cold /*K*/, double aspect_silo /*-*/, int silo_pairs /*-*/, double u_silo /*W/m^2-K*/,
		double & vol_one_temp_avail /*m3*/, double & vol_one_temp_total /*m3*/, double & d_silo /*m*/, double& h_silo /*m*/,
		double & q_dot_loss_des /*MWt*/  );


#endif   
