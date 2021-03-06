function spxy_tep_hessx = dess_2comp_echo1_exchg0_mag1_freq0_hessx(M0,ff,T1f,T1s,T2f,T2s,kap,flip,TR,TEp,R2pf,R2ps)
%DESS_2COMP_ECHO1_EXCHG0_MAG1_FREQ0_HESSX
%    SPXY_TEP_HESSX = DESS_2COMP_ECHO1_EXCHG0_MAG1_FREQ0_HESSX(M0,FF,T1F,T1S,T2F,T2S,KAP,FLIP,TR,TEP,R2PF,R2PS)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    22-Jan-2018 02:22:35
n = length(M0);
t0 = zeros(n,1);

t2 = 1.0./T1f;
t8 = TR.*t2;
t3 = exp(-t8);
t4 = flip.*kap;
t5 = cos(t4);
t6 = 1.0./T2f;
t25 = TR.*t6.*2.0;
t7 = exp(-t25);
t9 = -t3+t5;
t10 = t3.*t5;
t11 = t10-1.0;
t12 = 1.0./T1s;
t16 = TR.*t12;
t13 = exp(-t16);
t14 = 1.0./T2s;
t44 = TR.*t14.*2.0;
t15 = exp(-t44);
t17 = t5-t13;
t18 = t5.*t13;
t19 = t18-1.0;
t20 = flip.*kap.*(1.0./2.0);
t21 = tan(t20);
t22 = R2pf.*T2f;
t23 = t22+1.0;
t77 = TEp.*t6.*t23;
t24 = exp(-t77);
t26 = 1.0./t11;
t27 = -t7+1.0;
t28 = sqrt(t27);
t29 = t9.^2;
t30 = 1.0./t11.^2;
t31 = t7.*t29.*t30;
t32 = t31-1.0;
t33 = 1.0./t32;
t34 = sqrt(-t33);
t75 = t9.*t26.*t28.*t34;
t35 = -t75+1.0;
t36 = T1f.*4.0;
t37 = t5.^2;
t38 = exp(t25);
t39 = T1f.*2.0;
t40 = T2f+t39;
t41 = R2ps.*T2s;
t42 = t41+1.0;
t81 = TEp.*t14.*t42;
t43 = exp(-t81);
t45 = 1.0./t19;
t46 = -t15+1.0;
t47 = sqrt(t46);
t48 = t17.^2;
t49 = 1.0./t19.^2;
t50 = t15.*t48.*t49;
t51 = t50-1.0;
t52 = 1.0./t51;
t53 = sqrt(-t52);
t79 = t17.*t45.*t47.*t53;
t54 = -t79+1.0;
t55 = t37-1.0;
t56 = T1s.*4.0;
t57 = exp(t44);
t58 = T1s.*2.0;
t59 = T2s+t58;
t60 = R2pf+t6;
t139 = TEp.*t60;
t61 = exp(-t139);
t62 = t7-1.0;
t63 = t33.*t62;
t64 = sqrt(t63);
t140 = t9.*t26.*t64;
t65 = -t140+1.0;
t66 = 1.0./T2f.^2;
t67 = ff-1.0;
t68 = R2ps+t14;
t147 = TEp.*t68;
t69 = exp(-t147);
t70 = t15-1.0;
t71 = t52.*t70;
t72 = sqrt(t71);
t148 = t17.*t45.*t72;
t73 = -t148+1.0;
t74 = 1.0./T2s.^2;
t76 = abs(t35);
t78 = t21.*t24.*t76;
t80 = abs(t54);
t82 = t78-t21.*t43.*t80;
t83 = 1.0./T1f.^2;
t84 = sign(t35);
t85 = 1.0./sqrt(t27);
t86 = t38-1.0;
t87 = 1.0./sqrt(-t33);
t88 = T2f.*3.0;
t89 = t36+t88;
t90 = TR.*t2.*t6.*t89;
t91 = exp(t90);
t92 = T2f+t36;
t93 = TR.*t2.*t6.*t92;
t94 = exp(t93);
t95 = t37.*t94;
t96 = TR.*t2.*t6.*t40.*2.0;
t97 = exp(t96);
t155 = t5.*t97.*2.0;
t98 = t91+t95-t155;
t99 = T1f+T2f;
t100 = TR.*t2.*t6.*t99.*2.0;
t101 = exp(t100);
t102 = TR.*t2.*2.0;
t103 = exp(t102);
t104 = t37.*t38;
t105 = exp(t8);
t106 = t5.*t105.*2.0;
t107 = TR.*t2.*t6.*t40;
t108 = exp(t107);
t156 = t37.*t103;
t157 = t5.*t108.*2.0;
t109 = t101+t104+t106-t156-t157-1.0;
t110 = 1.0./t109.^2;
t111 = 1.0./T1s.^2;
t112 = sign(t54);
t113 = 1.0./sqrt(t46);
t114 = t57-1.0;
t115 = 1.0./sqrt(-t52);
t116 = T2s.*3.0;
t117 = t56+t116;
t118 = TR.*t12.*t14.*t117;
t119 = exp(t118);
t120 = T2s+t56;
t121 = TR.*t12.*t14.*t120;
t122 = exp(t121);
t123 = t37.*t122;
t124 = TR.*t12.*t14.*t59.*2.0;
t125 = exp(t124);
t214 = t5.*t125.*2.0;
t126 = t119+t123-t214;
t127 = T1s+T2s;
t128 = TR.*t12.*t14.*t127.*2.0;
t129 = exp(t128);
t130 = TR.*t12.*2.0;
t131 = exp(t130);
t132 = t37.*t57;
t133 = exp(t16);
t134 = t5.*t133.*2.0;
t135 = TR.*t12.*t14.*t59;
t136 = exp(t135);
t215 = t37.*t131;
t216 = t5.*t136.*2.0;
t137 = t129+t132+t134-t215-t216-1.0;
t138 = 1.0./t137.^2;
t141 = abs(t65);
t142 = sign(t65);
t143 = 1.0./sqrt(t63);
t144 = TR.*t7.*t33.*t66.*2.0;
t145 = 1.0./t32.^2;
t271 = TR.*t7.*t29.*t30.*t62.*t66.*t145.*2.0;
t146 = t144-t271;
t149 = abs(t73);
t150 = sign(t73);
t151 = 1.0./sqrt(t71);
t152 = TR.*t15.*t52.*t74.*2.0;
t153 = 1.0./t51.^2;
t335 = TR.*t15.*t48.*t49.*t70.*t74.*t153.*2.0;
t154 = t152-t335;
t158 = TR.*ff.*t7.*t21.*t24.*t55.*t83.*t84.*t85.*t86.*t87.*t98.*t110;
t159 = M0.*TR.*t7.*t21.*t24.*t55.*t83.*t84.*t85.*t86.*t87.*t98.*t110;
t160 = 1.0./T1f.^4;
t161 = TR.^2;
t162 = exp(-t102);
t163 = 1.0./T1f.^3;
t164 = t5-t105;
t165 = t5.*t105;
t166 = t165-1.0;
t167 = 1.0./t164.^3;
t168 = t55.^2;
t169 = 1.0./t109.^4;
t170 = T1f.*TEp;
t171 = T2f.*TR.*4.0;
t172 = R2pf.*T1f.*T2f.*TEp;
t173 = T1f.*TR.*6.0;
t174 = T2f.*TR.*6.0;
t175 = T1f.*TR.*4.0;
t176 = t170+t171+t172+t175;
t177 = t2.*t6.*t176;
t178 = exp(t177);
t179 = t170+t171+t172+t173;
t180 = t2.*t6.*t179;
t181 = exp(t180);
t182 = t170+t172+t173+t174;
t183 = t2.*t6.*t182;
t184 = exp(t183);
t185 = T2f.*TR.*3.0;
t186 = T2f.*TR.*5.0;
t187 = T1f.*TR.*8.0;
t188 = t170+t172+t175+t185;
t189 = t2.*t6.*t188;
t190 = exp(t189);
t191 = t170+t172+t175+t186;
t192 = t2.*t6.*t191;
t193 = exp(t192);
t194 = t170+t172+t173+t185;
t195 = t2.*t6.*t194;
t196 = exp(t195);
t197 = t170+t172+t173+t186;
t198 = t2.*t6.*t197;
t199 = exp(t198);
t200 = T2f.*TR.*2.0;
t201 = t37.^2;
t202 = t170+t172+t173+t200;
t203 = t2.*t6.*t202;
t204 = exp(t203);
t205 = t170+t172+t175+t200;
t206 = t2.*t6.*t205;
t207 = exp(t206);
t208 = t170+t172+t174+t175;
t209 = t2.*t6.*t208;
t210 = exp(t209);
t211 = dirac(t35);
t212 = t164.^2;
t213 = 1.0./t62;
t217 = 1.0./T1s.^4;
t218 = exp(-t130);
t219 = 1.0./T1s.^3;
t220 = TR.*t5;
t221 = t5-t133;
t222 = t5.*t133;
t223 = t222-1.0;
t224 = 1.0./t221.^3;
t225 = TR.*4.0;
t226 = 1.0./t137.^4;
t227 = T1s.*TEp;
t228 = T2s.*TR.*4.0;
t229 = R2ps.*T1s.*T2s.*TEp;
t230 = T1s.*TR.*6.0;
t231 = T2s.*TR.*6.0;
t232 = T1s.*TR.*4.0;
t233 = t227+t228+t229+t232;
t234 = t12.*t14.*t233;
t235 = exp(t234);
t236 = t227+t228+t229+t230;
t237 = t12.*t14.*t236;
t238 = exp(t237);
t239 = t227+t229+t230+t231;
t240 = t12.*t14.*t239;
t241 = exp(t240);
t242 = T2s.*TR.*3.0;
t243 = T2s.*TR.*5.0;
t244 = T1s.*TR.*8.0;
t245 = t227+t229+t232+t242;
t246 = t12.*t14.*t245;
t247 = exp(t246);
t248 = t227+t229+t232+t243;
t249 = t12.*t14.*t248;
t250 = exp(t249);
t251 = t227+t229+t230+t242;
t252 = t12.*t14.*t251;
t253 = exp(t252);
t254 = t227+t229+t230+t243;
t255 = t12.*t14.*t254;
t256 = exp(t255);
t257 = T2s.*TR.*2.0;
t258 = t227+t229+t230+t257;
t259 = t12.*t14.*t258;
t260 = exp(t259);
t261 = t227+t229+t232+t257;
t262 = t12.*t14.*t261;
t263 = exp(t262);
t264 = t227+t229+t231+t232;
t265 = t12.*t14.*t264;
t266 = exp(t265);
t267 = dirac(t54);
t268 = t221.^2;
t269 = 1.0./t70;
t270 = TEp.*ff.*t21.*t61.*t66.*t141;
t272 = t270-ff.*t9.*t21.*t26.*t61.*t142.*t143.*t146.*(1.0./2.0);
t273 = M0.*TEp.*t21.*t61.*t66.*t141;
t274 = t273-M0.*t9.*t21.*t26.*t61.*t142.*t143.*t146.*(1.0./2.0);
t275 = T1f.*TEp.*2.0;
t276 = T1f.*TR.*2.0;
t277 = T2f.*TR;
t278 = R2pf.*T1f.*T2f.*TEp.*2.0;
t279 = t275+t276+t277+t278;
t280 = exp(-t2.*t6.*t279);
t281 = 1.0./t27.^(3.0./2.0);
t282 = TEp.*t181;
t283 = TEp.*t184;
t284 = t170+t172+t174+t187;
t285 = t2.*t6.*t284;
t286 = exp(t285);
t287 = TR.*t178.*2.0;
t288 = TR.*t184;
t289 = TEp.*t5.*t190.*2.0;
t290 = TEp.*t5.*t193.*2.0;
t291 = t170+t172+t186+t187;
t292 = t2.*t6.*t291;
t293 = exp(t292);
t294 = TEp.*t5.*t293.*4.0;
t295 = TR.*t5.*t196.*6.0;
t296 = TR.*t5.*t199.*2.0;
t297 = TEp.*t5.*t37.*t190.*2.0;
t298 = TEp.*t37.*t204;
t299 = TEp.*t201.*t204;
t300 = TEp.*t5.*t37.*t193.*2.0;
t301 = TEp.*t37.*t181.*1.0e1;
t302 = TEp.*t181.*t201;
t303 = t170+t172+t187+t200;
t304 = t2.*t6.*t303;
t305 = exp(t304);
t306 = t170+t172+t185+t187;
t307 = t2.*t6.*t306;
t308 = exp(t307);
t309 = TEp.*t5.*t37.*t308.*4.0;
t310 = TEp.*t37.*t184;
t311 = t170+t171+t172+t187;
t312 = t2.*t6.*t311;
t313 = exp(t312);
t314 = TR.*t37.*t207.*2.0;
t315 = TR.*t37.*t178.*8.0;
t316 = TR.*t178.*t201.*2.0;
t317 = TR.*t201.*t204;
t318 = TR.*t5.*t37.*t196.*2.0;
t319 = TR.*t37.*t210.*2.0;
t320 = TR.*t5.*t37.*t199.*6.0;
t321 = t282+t283+t287+t288+t289+t290+t294+t295+t296+t297+t298+t299+t300+t301+t302+t309+t310+t314+t315+t316+t317+t318+t319+t320-TEp.*t178-TEp.*t286-TR.*t181.*3.0-TEp.*t5.*t196.*2.0-TEp.*t5.*t199.*6.0-TEp.*t37.*t178.*4.0-TEp.*t37.*t207-TEp.*t37.*t210-TEp.*t37.*t313.*6.0-TEp.*t178.*t201-TEp.*t201.*t305-TR.*t5.*t190.*4.0-TR.*t5.*t193.*4.0-TR.*t37.*t181.*6.0-TR.*t37.*t184.*3.0-TR.*t37.*t204.*3.0-TR.*t181.*t201.*3.0-TEp.*t5.*t37.*t196.*6.0-TEp.*t5.*t37.*t199.*2.0-TR.*t5.*t37.*t190.*4.0-TR.*t5.*t37.*t193.*4.0;
t322 = M0.*TR.*ff.*t21.*t32.*t55.*t66.*t83.*t84.*t86.*t87.*t169.*t212.*t280.*t281.*t321;
t323 = t103-1.0;
t324 = t322-M0.*ff.*t9.*t21.*t24.*t26.*t32.*t66.*t83.*t86.*t98.*t161.*t168.*t169.*t211.*t212.*t213.*t323.*2.0;
t325 = 1.0./t84;
t326 = t75-1.0;
t327 = 1.0./T2f.^4;
t328 = R2pf.*T2f.*TEp;
t329 = TEp-t225+t328;
t330 = exp(-t6.*t329);
t331 = t212.^2;
t332 = t323.^2;
t395 = TR.*2.0;
t333 = TEp+t328-t395;
t334 = exp(-t6.*t333);
t336 = t17.*t21.*t45.*t67.*t69.*t150.*t151.*t154.*(1.0./2.0);
t337 = t336-TEp.*t21.*t67.*t69.*t74.*t149;
t338 = M0.*t17.*t21.*t45.*t69.*t150.*t151.*t154.*(1.0./2.0);
t339 = t338-M0.*TEp.*t21.*t69.*t74.*t149;
t340 = T1s.*TEp.*2.0;
t341 = T1s.*TR.*2.0;
t342 = T2s.*TR;
t343 = R2ps.*T1s.*T2s.*TEp.*2.0;
t344 = t340+t341+t342+t343;
t345 = exp(-t12.*t14.*t344);
t346 = 1.0./t46.^(3.0./2.0);
t347 = TEp.*t238;
t348 = TEp.*t241;
t349 = t227+t229+t231+t244;
t350 = t12.*t14.*t349;
t351 = exp(t350);
t352 = TR.*t235.*2.0;
t353 = TR.*t241;
t354 = TEp.*t5.*t247.*2.0;
t355 = TEp.*t5.*t250.*2.0;
t356 = t227+t229+t243+t244;
t357 = t12.*t14.*t356;
t358 = exp(t357);
t359 = TEp.*t5.*t358.*4.0;
t360 = TR.*t5.*t253.*6.0;
t361 = TR.*t5.*t256.*2.0;
t362 = TEp.*t5.*t37.*t247.*2.0;
t363 = TEp.*t37.*t260;
t364 = TEp.*t201.*t260;
t365 = TEp.*t5.*t37.*t250.*2.0;
t366 = TEp.*t37.*t238.*1.0e1;
t367 = TEp.*t201.*t238;
t368 = t227+t229+t244+t257;
t369 = t12.*t14.*t368;
t370 = exp(t369);
t371 = t227+t229+t242+t244;
t372 = t12.*t14.*t371;
t373 = exp(t372);
t374 = TEp.*t5.*t37.*t373.*4.0;
t375 = TEp.*t37.*t241;
t376 = t227+t228+t229+t244;
t377 = t12.*t14.*t376;
t378 = exp(t377);
t379 = TR.*t37.*t263.*2.0;
t380 = TR.*t37.*t235.*8.0;
t381 = TR.*t201.*t235.*2.0;
t382 = TR.*t201.*t260;
t383 = TR.*t5.*t37.*t253.*2.0;
t384 = TR.*t37.*t266.*2.0;
t385 = TR.*t5.*t37.*t256.*6.0;
t386 = t347+t348+t352+t353+t354+t355+t359+t360+t361+t362+t363+t364+t365+t366+t367+t374+t375+t379+t380+t381+t382+t383+t384+t385-TEp.*t235-TEp.*t351-TR.*t238.*3.0-TEp.*t5.*t253.*2.0-TEp.*t5.*t256.*6.0-TEp.*t37.*t235.*4.0-TEp.*t37.*t263-TEp.*t37.*t266-TEp.*t37.*t378.*6.0-TEp.*t201.*t235-TEp.*t201.*t370-TR.*t5.*t247.*4.0-TR.*t5.*t250.*4.0-TR.*t37.*t238.*6.0-TR.*t37.*t241.*3.0-TR.*t37.*t260.*3.0-TR.*t201.*t238.*3.0-TEp.*t5.*t37.*t253.*6.0-TEp.*t5.*t37.*t256.*2.0-TR.*t5.*t37.*t247.*4.0-TR.*t5.*t37.*t250.*4.0;
t387 = t131-1.0;
t388 = M0.*t17.*t21.*t43.*t45.*t51.*t67.*t74.*t111.*t114.*t126.*t161.*t168.*t226.*t267.*t268.*t269.*t387.*2.0;
t389 = t388-M0.*TR.*t21.*t51.*t55.*t67.*t74.*t111.*t112.*t114.*t115.*t226.*t268.*t345.*t346.*t386;
t390 = TEp.^2;
t391 = 1.0./t112;
t392 = t79-1.0;
t393 = 1.0./T2s.^4;
t394 = R2ps.*T2s.*TEp;
t396 = TEp-t225+t394;
t397 = exp(-t14.*t396);
t398 = t268.^2;
t399 = t387.^2;
t400 = TEp+t394-t395;
t401 = exp(-t14.*t400);
% spxy_tep_hessx = reshape([0.0,t82,t158,-TR.*t15.*t21.*t43.*t55.*t67.*t111.*t112.*t113.*t114.*t115.*t126.*t138,t272,t337,t82,0.0,t159,-M0.*TR.*t15.*t21.*t43.*t55.*t111.*t112.*t113.*t114.*t115.*t126.*t138,t274,t339,t158,t159,-M0.*ff.*t21.*t24.*t84.*(-t3.*t26.*t28.*t34.*t160.*t161+TR.*t3.*t26.*t28.*t34.*t163.*2.0+t5.*t28.*t30.*t34.*t160.*t161.*t162.*2.0+t9.*1.0./t11.^3.*t28.*t34.*t37.*t160.*t161.*t162.*2.0+TR.*t3.*t5.*t9.*t28.*t30.*t34.*t163.*2.0-t3.*t5.*t9.*t28.*t30.*t34.*t160.*t161+t7.*t26.*t28.*t55.*t87.*t145.*t160.*t161.*t166.*t167.*2.0-t9.*t26.*t28.*1.0./t32.^3.*t87.*t160.*t161.*1.0./t164.^6.*t166.^2.*t168.*exp(TR.*t2.*t6.*(T2f-t39).*2.0).*3.0-(TR.*t9.*t26.*t28.*t55.*t87.*t105.*t145.*t160.*(t220+T1f.*t5.*2.0-T1f.*t105.*2.0+TR.*t105.*2.0+T1f.*t5.*t103.*2.0-T1f.*t37.*t105.*2.0-TR.*t5.*t103-TR.*t37.*t105.*2.0))./(exp(TR.*t2.*t6.*(T1f+T2f.*2.0).*2.0)-t5.*exp(TR.*t2.*t6.*(t39+t88)).*4.0+t37.*t101.*6.0+t38.*t201-t5.*t37.*t108.*4.0)+t5.*t7.*t9.*t28.*t30.*t55.*t87.*t145.*t160.*t161.*t166.*t167.*2.0)+M0.*ff.*t21.*t32.*t86.^2.*t98.^2.*t160.*t161.*t168.*t169.*t211.*t213.*exp(-t6.*(TEp+t225+t328)).*2.0,0.0,t324,0.0,-TR.*t15.*t21.*t43.*t55.*t67.*t111.*t112.*t113.*t114.*t115.*t126.*t138,-M0.*TR.*t15.*t21.*t43.*t55.*t111.*t112.*t113.*t114.*t115.*t126.*t138,0.0,M0.*t21.*t43.*t67.*t112.*(-t13.*t45.*t47.*t53.*t161.*t217+TR.*t13.*t45.*t47.*t53.*t219.*2.0+t5.*t47.*t49.*t53.*t161.*t217.*t218.*2.0+t17.*1.0./t19.^3.*t37.*t47.*t53.*t161.*t217.*t218.*2.0+TR.*t5.*t13.*t17.*t47.*t49.*t53.*t219.*2.0-t5.*t13.*t17.*t47.*t49.*t53.*t161.*t217+t15.*t45.*t47.*t55.*t115.*t153.*t161.*t217.*t223.*t224.*2.0-t17.*t45.*t47.*1.0./t51.^3.*t115.*t161.*t168.*t217.*1.0./t221.^6.*t223.^2.*exp(TR.*t12.*t14.*(T2s-t58).*2.0).*3.0-(TR.*t17.*t45.*t47.*t55.*t115.*t133.*t153.*t217.*(t220+T1s.*t5.*2.0-T1s.*t133.*2.0+TR.*t133.*2.0+T1s.*t5.*t131.*2.0-T1s.*t37.*t133.*2.0-TR.*t5.*t131-TR.*t37.*t133.*2.0))./(exp(TR.*t12.*t14.*(T1s+T2s.*2.0).*2.0)-t5.*exp(TR.*t12.*t14.*(t58+t116)).*4.0+t37.*t129.*6.0+t57.*t201-t5.*t37.*t136.*4.0)+t5.*t15.*t17.*t47.*t49.*t55.*t115.*t153.*t161.*t217.*t223.*t224.*2.0)-M0.*t21.*t51.*t67.*t114.^2.*t126.^2.*t161.*t168.*t217.*t226.*t267.*t269.*exp(-t14.*(TEp+t225+t394)).*2.0,0.0,t389,t272,t274,t324,0.0,-M0.*ff.*t21.*t24.*t325.*t326.*t327.*t390+M0.*1.0./T2f.^3.*TEp.*ff.*t21.*t24.*t325.*t326.*2.0+M0.*ff.*t21.*t29.*t30.*t32.*t161.*t168.*t169.*t211.*t213.*t327.*t330.*t331.*t332.*2.0-M0.*TEp.*TR.*ff.*t9.*t21.*t26.*t55.*t84.*t85.*t87.*t110.*t212.*t323.*t327.*t334.*2.0-M0.*TR.*ff.*t9.*t21.*t26.*t55.*t84.*t85.*t87.*1.0./t109.^3.*t212.*t323.*t327.*t334.*(T2f+TR-T2f.*t101+TR.*t101-T2f.*t37.*t38-T2f.*t5.*t105.*2.0+T2f.*t5.*t108.*2.0+T2f.*t37.*t103+TR.*t37.*t38-TR.*t5.*t105.*2.0-TR.*t5.*t108.*2.0+TR.*t37.*t103).*2.0-M0.*ff.*t9.*t21.*t26.*t32.*t84.*t87.*t161.*t168.*t169.*t281.*t327.*t330.*t331.*t332,0.0,t337,t339,0.0,t389,0.0,M0.*t21.*t43.*t67.*t390.*t391.*t392.*t393-M0.*1.0./T2s.^3.*TEp.*t21.*t43.*t67.*t391.*t392.*2.0-M0.*t21.*t48.*t49.*t51.*t67.*t161.*t168.*t226.*t267.*t269.*t393.*t397.*t398.*t399.*2.0+M0.*TEp.*TR.*t17.*t21.*t45.*t55.*t67.*t112.*t113.*t115.*t138.*t268.*t387.*t393.*t401.*2.0+M0.*TR.*t17.*t21.*t45.*t55.*t67.*t112.*t113.*t115.*1.0./t137.^3.*t268.*t387.*t393.*t401.*(T2s+TR-T2s.*t129+TR.*t129-T2s.*t37.*t57-T2s.*t5.*t133.*2.0+T2s.*t5.*t136.*2.0+T2s.*t37.*t131+TR.*t37.*t57-TR.*t5.*t133.*2.0-TR.*t5.*t136.*2.0+TR.*t37.*t131).*2.0+M0.*t17.*t21.*t45.*t51.*t67.*t112.*t115.*t161.*t168.*t226.*t346.*t393.*t397.*t398.*t399],[1,6,6]);
spxy_tep_hessx = reshape([t0,t82,t158,-TR.*t15.*t21.*t43.*t55.*t67.*t111.*t112.*t113.*t114.*t115.*t126.*t138,t272,t337,t82,t0,t159,-M0.*TR.*t15.*t21.*t43.*t55.*t111.*t112.*t113.*t114.*t115.*t126.*t138,t274,t339,t158,t159,-M0.*ff.*t21.*t24.*t84.*(-t3.*t26.*t28.*t34.*t160.*t161+TR.*t3.*t26.*t28.*t34.*t163.*2.0+t5.*t28.*t30.*t34.*t160.*t161.*t162.*2.0+t9.*1.0./t11.^3.*t28.*t34.*t37.*t160.*t161.*t162.*2.0+TR.*t3.*t5.*t9.*t28.*t30.*t34.*t163.*2.0-t3.*t5.*t9.*t28.*t30.*t34.*t160.*t161+t7.*t26.*t28.*t55.*t87.*t145.*t160.*t161.*t166.*t167.*2.0-t9.*t26.*t28.*1.0./t32.^3.*t87.*t160.*t161.*1.0./t164.^6.*t166.^2.*t168.*exp(TR.*t2.*t6.*(T2f-t39).*2.0).*3.0-(TR.*t9.*t26.*t28.*t55.*t87.*t105.*t145.*t160.*(t220+T1f.*t5.*2.0-T1f.*t105.*2.0+TR.*t105.*2.0+T1f.*t5.*t103.*2.0-T1f.*t37.*t105.*2.0-TR.*t5.*t103-TR.*t37.*t105.*2.0))./(exp(TR.*t2.*t6.*(T1f+T2f.*2.0).*2.0)-t5.*exp(TR.*t2.*t6.*(t39+t88)).*4.0+t37.*t101.*6.0+t38.*t201-t5.*t37.*t108.*4.0)+t5.*t7.*t9.*t28.*t30.*t55.*t87.*t145.*t160.*t161.*t166.*t167.*2.0)+M0.*ff.*t21.*t32.*t86.^2.*t98.^2.*t160.*t161.*t168.*t169.*t211.*t213.*exp(-t6.*(TEp+t225+t328)).*2.0,t0,t324,t0,-TR.*t15.*t21.*t43.*t55.*t67.*t111.*t112.*t113.*t114.*t115.*t126.*t138,-M0.*TR.*t15.*t21.*t43.*t55.*t111.*t112.*t113.*t114.*t115.*t126.*t138,t0,M0.*t21.*t43.*t67.*t112.*(-t13.*t45.*t47.*t53.*t161.*t217+TR.*t13.*t45.*t47.*t53.*t219.*2.0+t5.*t47.*t49.*t53.*t161.*t217.*t218.*2.0+t17.*1.0./t19.^3.*t37.*t47.*t53.*t161.*t217.*t218.*2.0+TR.*t5.*t13.*t17.*t47.*t49.*t53.*t219.*2.0-t5.*t13.*t17.*t47.*t49.*t53.*t161.*t217+t15.*t45.*t47.*t55.*t115.*t153.*t161.*t217.*t223.*t224.*2.0-t17.*t45.*t47.*1.0./t51.^3.*t115.*t161.*t168.*t217.*1.0./t221.^6.*t223.^2.*exp(TR.*t12.*t14.*(T2s-t58).*2.0).*3.0-(TR.*t17.*t45.*t47.*t55.*t115.*t133.*t153.*t217.*(t220+T1s.*t5.*2.0-T1s.*t133.*2.0+TR.*t133.*2.0+T1s.*t5.*t131.*2.0-T1s.*t37.*t133.*2.0-TR.*t5.*t131-TR.*t37.*t133.*2.0))./(exp(TR.*t12.*t14.*(T1s+T2s.*2.0).*2.0)-t5.*exp(TR.*t12.*t14.*(t58+t116)).*4.0+t37.*t129.*6.0+t57.*t201-t5.*t37.*t136.*4.0)+t5.*t15.*t17.*t47.*t49.*t55.*t115.*t153.*t161.*t217.*t223.*t224.*2.0)-M0.*t21.*t51.*t67.*t114.^2.*t126.^2.*t161.*t168.*t217.*t226.*t267.*t269.*exp(-t14.*(TEp+t225+t394)).*2.0,t0,t389,t272,t274,t324,t0,-M0.*ff.*t21.*t24.*t325.*t326.*t327.*t390+M0.*1.0./T2f.^3.*TEp.*ff.*t21.*t24.*t325.*t326.*2.0+M0.*ff.*t21.*t29.*t30.*t32.*t161.*t168.*t169.*t211.*t213.*t327.*t330.*t331.*t332.*2.0-M0.*TEp.*TR.*ff.*t9.*t21.*t26.*t55.*t84.*t85.*t87.*t110.*t212.*t323.*t327.*t334.*2.0-M0.*TR.*ff.*t9.*t21.*t26.*t55.*t84.*t85.*t87.*1.0./t109.^3.*t212.*t323.*t327.*t334.*(T2f+TR-T2f.*t101+TR.*t101-T2f.*t37.*t38-T2f.*t5.*t105.*2.0+T2f.*t5.*t108.*2.0+T2f.*t37.*t103+TR.*t37.*t38-TR.*t5.*t105.*2.0-TR.*t5.*t108.*2.0+TR.*t37.*t103).*2.0-M0.*ff.*t9.*t21.*t26.*t32.*t84.*t87.*t161.*t168.*t169.*t281.*t327.*t330.*t331.*t332,t0,t337,t339,t0,t389,t0,M0.*t21.*t43.*t67.*t390.*t391.*t392.*t393-M0.*1.0./T2s.^3.*TEp.*t21.*t43.*t67.*t391.*t392.*2.0-M0.*t21.*t48.*t49.*t51.*t67.*t161.*t168.*t226.*t267.*t269.*t393.*t397.*t398.*t399.*2.0+M0.*TEp.*TR.*t17.*t21.*t45.*t55.*t67.*t112.*t113.*t115.*t138.*t268.*t387.*t393.*t401.*2.0+M0.*TR.*t17.*t21.*t45.*t55.*t67.*t112.*t113.*t115.*1.0./t137.^3.*t268.*t387.*t393.*t401.*(T2s+TR-T2s.*t129+TR.*t129-T2s.*t37.*t57-T2s.*t5.*t133.*2.0+T2s.*t5.*t136.*2.0+T2s.*t37.*t131+TR.*t37.*t57-TR.*t5.*t133.*2.0-TR.*t5.*t136.*2.0+TR.*t37.*t131).*2.0+M0.*t17.*t21.*t45.*t51.*t67.*t112.*t115.*t161.*t168.*t226.*t346.*t393.*t397.*t398.*t399],[n,6,6]);
