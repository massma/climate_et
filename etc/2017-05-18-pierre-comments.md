### Pierre ###

I quickly played to see how ET will change in the future.
I modified Rnet, esat(T) and VPD and the stomatal response in Penman Monteith.
So even for the full ET it is all VPD driven…
I am surprised because esat and VPD have the same exponential dependence

I assumed RCP8.5 warming (3.7K for up to 1200 ppm of CO2 equivalent). Then I can plot the réponse of radiation vs temperature as a function of CO2

So now I am pretty sure that the increase in ET has nothing to do with VPD changes (opposite to what Jack is saying below)
This is interesting

Here is the reason:

ET =( Delta (Rn-G) + RCP/ra*VPD ) / (Delta + gamma*(1+rs/ra)) based on penman monteith

divide by Delta the num and den

ET =( (Rn-G) + RCP/ra*VPD/Delta ) / (1 + gamma/Delta*(1+rs/ra)) based on penman monteith

VPD/Delta does not change (same exponential dependence of both d Esat/dt = Esat * stuff)
Rn changes by 8.5 W/m^2 in RCP 8.5

The main dependence is the gamma/Delta term: in other words the surface becomes less and less important with global warming. Makes sense as we are decreasing the radiation dependence on the surface vapor pressure deficit across the leaf.

Definitely worth writing up.

### Response ###

Still some parts I do not understand (sorry if I am slow).... Why would VPD/delta scale as 1 with increasing T? If VPD scales as e_sat(T), but delta scales as d(e_sat(T))/dT, then wouldn't they scale differently with increasing temperature? Also, if VPD doesn't matter how come the idealized penman monteith plots you made show that most of the change in ET is due to VPD? Thanks for your help getting me to understand this.

### Pierre's clarification ###

In fact desat/dT is proportional to Esat times a very minorly T dependent function. 
Indeed the VPD plot appears to be strong but in this one I didn't include the T dependence in desat/dT. If you include it it cancels out the VPD term and then only appears as 1/Delta rs/ra.

