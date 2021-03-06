TITLE fast activated potassium Kv3 (Kv3.1/3.4) channel for GPe neuron

COMMENT
 modeled by Gunay et al., 2008
 implemented in NEURON by Kitano, 2011
ENDCOMMENT

UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
}

NEURON {
 SUFFIX Kv3
 USEION k READ ek WRITE ik
 RANGE gkbar, gk, ik
 GLOBAL scale_a
}

PARAMETER {
 v (mV)
 dt (ms)
 gkbar  = 0.001 (S/cm2)
 ek (mV)

 theta_m = -26.0 (mV)
 k_m = 7.8 (mV)
 tau_m0 = 0.1 (ms)
 tau_m1 = 14.0 (ms)
 phi_m = -26.0 (mV)
 sigma_m0 = 13.0 (mV)
 sigma_m1 = -12.0 (mV)

 h0 = 0.6
 theta_h = -20.0 (mV)
 k_h = -10.0 (mV)
 tau_h0 = 7.0 (ms)
 tau_h1 = 33.0 (ms)
 phi_h = 0.0 (mV)
 sigma_h0 = 10.0 (mV)
 sigma_h1 = -10.0 (mV)
 
 scale_a = 1 (1)
}

STATE {
 m h
}

ASSIGNED { 
 ik (mA/cm2)
 minf
 taum (ms)
 hinf
 tauh (ms)
 gk (S/cm2)
}

BREAKPOINT {
 SOLVE states METHOD cnexp
 gk  = gkbar*m*m*m*m*h
 ik = gk * (v-ek)
}

UNITSOFF

INITIAL {
 settables(v)
 m = minf
 h = hinf
}

DERIVATIVE states {  
 settables(v)
 m' = (minf - m)/taum
 h' = (hinf - h)/tauh
}

PROCEDURE settables(v) {
        TABLE minf, taum, hinf, tauh FROM -100 TO 100 WITH 400

	minf = 1.0 / (1.0 + exp((theta_m - v)/k_m))
	taum = tau_m0 + (tau_m1 - tau_m0)/(scale_a*(exp((phi_m - v)/sigma_m0) + exp((phi_m - v)/sigma_m1)))
	hinf = h0 + (1.0 - h0) / (1.0 + exp((theta_h - v)/k_h))
	tauh = tau_h0 + (tau_h1 - tau_h0)/(exp((phi_h - v)/sigma_h0) + exp((phi_h - v)/sigma_h1))
}

UNITSON






