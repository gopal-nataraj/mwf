function ssxy_te = spgr_2comp_exchg0_mag1_freq0(M0,ff,T1f,T1s,T2f,T2s,kap,flip,TR,TE,R2pf,R2ps)
%SPGR_2COMP_EXCHG0_MAG1_FREQ0
%    SSXY_TE = SPGR_2COMP_EXCHG0_MAG1_FREQ0(M0,FF,T1F,T1S,T2F,T2S,KAP,FLIP,TR,TE,R2PF,R2PS)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    21-Jan-2018 18:21:48

t2 = flip.*kap;
t3 = 1.0./T1s;
t4 = TR.*t3;
t5 = exp(t4);
t6 = sin(t2);
t7 = cos(t2);
t8 = 1.0./T1f;
t9 = TR.*t8;
t10 = exp(t9);
ssxy_te = (M0.*ff.*t6.*t10.*exp(-TE.*(R2pf+1.0./T2f)).*(exp(-t9)-1.0))./(t7-t10)+(M0.*t5.*t6.*exp(-TE.*(R2ps+1.0./T2s)).*(exp(-t4)-1.0).*(ff-1.0))./(t5-t7);
