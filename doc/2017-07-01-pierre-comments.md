### Pierre: ###

I think in the summer you could start writing a paper on temperature and VPD influence on evaporation (like what we discussed). I think everyone is confused as I mentioned and it is really easy to now write a paper because you have Changjie’s results. I know someone who might want to write something on VPD and evap so it might be good to rush this one (while you are progressing on the Amazon).

### Response: ###

I was just reviewing Jack Scheff's papers, and I'm not sure how what we are proposing is very different from this one:

http://www.ldeo.columbia.edu/~jscheff/SF2013ScalingPotEvap.pdf

They also conclude that changes are mostly due to temperature effect... how would our plan be very different?

### Pierre's clarification ###

I had forgotten that Jack went so well into the theory and explanation on the CC scaling and Rn-G. I really like most of his work. Anyhow, potential evaporation is to me more and more a useless concept. In reality you will never observe the type of dependence on CC scaling he is highlighting. The reason for this is that you also have a VPD dependence in stomata conductance. This is what I meant: if you disconnect this g_sto term in ET indeed you will observe a strong T dependence as it impacts both VPD and esat. In reality the dependence will be very mild and only due to the interaction of the slope with the g_sto VPD dependence. In other words, we will never observe the rise jack talks about with rising T. It will be very mild. It also means that people exaggerate the importance of VPD for ecosystem functioning and ET in particular.

I think we can nicely link this to Alexis Berg’s recent paper on the two levels of soil smoiture. He also has a paper on increased aridity, but again to me it is pretty useless as actual ET will never be working this way. BUT since Jack explained many of the terms well we need to explain the ET story in a nice and consistent way. I think a nice paper could be that we could explain the surface drying compared to near steady deeper soil moisture by this type of analysis.
Surface soil moisture is pretty much at potential evaporation rate (times a soil moisture stress). Deeper rooting depth soil moisture is transpiration dominated.  The two are going to evolve very differently in the future. The surface one indeed exhibits a CC scaling. The bottom one nearly has no dependence and only the CO2 effect is important (you can use the atmospheric warming experiment only to see this). So that could be a (nice) first paper.

The second paper I had in mind is whether VPD is a driver or reducer of ET and GPP. To do this use Penman Monteith (Matlab script I sent) and investigate the impact of the numerator vs denominator (g_sto dependence). The challenge is that you need a (good) model for photosynthesis and g_st as g_c = LAI.g_st, with g_st=gst_min + m GPP/VPD^n. The n dependence will be critical (except if it cancels out see below), with n=0.5 (Medlyn optimal model or 1 for Luening’s model, we should try both). m is a water use efficiency, ecosystem dependent. Now the challenge is that GPP is correlated with ET. I thiink there is a trick. uWUE (underlying water use efficiency) has been shown to be quite conserved across time scales (I will send you all papers I am talking about in a second email) uWUE = GPP.VPD^n/ET. If you use n=1 then you need to use intrinsic water use efficiency also relatively well conserved. Then plug that into PM and solve for ET again. Once you have ET you also have GPP (I would guess analytically). Then you can solve the air temperature and VPD dependence. To check only the VPD dependence just change relative humidity but not T. Then plot the influence of either the demand term (numerator in PM) or VPD dependence at the denominator. Then you can conclude.


