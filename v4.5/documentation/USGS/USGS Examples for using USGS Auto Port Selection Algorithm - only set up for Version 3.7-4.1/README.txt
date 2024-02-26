README.txt                                                     11-December-2014

Example Applications with Modified version 3.7/4.2 blending algorithm by USGS

This package contains four example applications of CE-QUAL-W2 models, meant to
illustrate the features of the modified blending algorithm.  

--------------------------------------------------------------------------------
Example #1:  det_normal_uro-float_400fmin

This is a model of Detroit Lake, OR, USA, for a "normal" water year, using a combination
of hypothetical outlets, where one of them is a floating outlet with a 400 cfs
minimum flow (float_400fmin).

This scenario uses 4 outlets: 
 1: floating weir, priority 1, 2.3 m depth, minimum 400 cfs, maximum 5600 cfs
 2: spillway, priority -1 (nonblended)
 3: lower power outlet, priority 1, maximum 5600 cfs
 4: regulating outlet, priority -1 (nonblended), maximum 5600 cfs

Outlets 2 and 4 represent outlets used to "spill" excess flow that the other
outlets cannot handle.  Those flows are preset in the outflow file, but the
temperature effects of those releases are accounted for by including them in the
blending group but giving them priorities of -1, which tells the model to
account for that heat but not adjust the flows through those two outlets.  The
model estimates the temperatures released by those flows through the priority -1
outlets, then adjusts the temperature target for the blended releases accordingly.
In this case, we have two outlets of the same priority balancing flows from a
fixed elevation (outlet #3) and from a floating outlet (#1).  It's not the most
exciting example, but it's a good example to illustrate the "nonblended" outflows
that are still accounted for in the temperature blending calculations.

The temperature target is in the file dynsplit_selective1.npt.  It was set up to
try to eject lots of heat through mid-summer, then concentrate on colder water
releases based on a no-dams assessment of temperatures later in the year.

Note that when the maximum flow of an outlet is exceeded, even for a non-blended
outlet (priority -1), the maximum flow criterion is honored and excess releases
are shifted to other outlets.

--------------------------------------------------------------------------------
Example #2:  lop-dex_lopFloat_20ppmin

This is a model for calendar year 2002 for Lookout Point and Dexter Lakes
(lop-dex) using a floating outlet at Lookout Point (branch 1, outlet 1, 1-m
depth) with a 20% minimum power production constraint at Lookout Point Dam
(20ppmin).

Branch 1 is Lookout Point Lake, and the dam is given 3 outlets:
 1: floating outlet at 1-m depth, priority 2, no minimum or maximum flow constraint
 2: power outlet, priority 1, 20% minimum flow constraint
 3: regulating outlet, priority 2, maximum head constraint of 51.42 m
TSSHARE is set to OFF, which tells the model to decide which of the two
priority-2 outlets to use in blending with the power outlet.  There are times
when the RO (#3) is not available because of the maximum head constraint.  In
that case, blending occurs between the power outlet and the floating outlet. 
Late in the year, when we need cold-water releases and the lake level is lower,
blending occurs between the power outlet and the RO. In this scenario, at least
20% of the releases at Lookout Point Dam are constrained to go through the power
outlet (minfrac=0.2).

Branch 2 is Dexter Lake.  Two outlets.  Nothing special.  Spillway is outlet #1
and has priority 2.  Power outlet is #2 and has priority 1.

Temperature targets are in dynsplit_selective1.npt for Lookout Point and
dynsplit_selective2.npt for Dexter.  They happen to be the same.

--------------------------------------------------------------------------------
Example #3:  det_multigate_example3

This example is a hypothetical case based on the Detroit Lake example, in which
a multiple-gate tower of 8 outlets is used to blend releases to meet a release
temperature target.  The outlets are all fixed-elevation and arranged from 476 m
to 385 m so that each is 13 meters apart, vertically.

The intention is to blend releases from the lowest and the highest available
outlets.  Therefore, priorities are set so that the deepest outlet has a priority
of 1, and the other outlets, from shallowest to deepest, have priorities from 2
to 8.  In this way, the deepest outlet will be blended with the one that is
nearest to the surface.  Each outlet is given a minimum head constraint of 2 m,
so that the outlet cannot be used if it is too near the water surface.  No other
flow or head constraints are specified.

Note that the results from this scenario would be similar to those resulting from
setting all of the outlet priorities to the same number.  In that case, the model
would first fulfill any minimum flow constraints (none here), then blend between
the lowest and highest available outlets.  That is pretty much what is done in
this example, only through the explicit setting of priorities.

--------------------------------------------------------------------------------
Example #4:  det_multigate_example4

This example is similar to example #3, except that some other constraints are
added and a different way of dealing with priorities is used.

In this example, outlet #6 is given a priority of 1 and all of the other outlets
are given a priority of 2.  TSSHARE is set to OFF, meaning that the model is
supposed to choose one of the priority-2 outlets to blend releases with the
single priority 1 outlet.  The choice is made after fulfilling any other minimum
flow constraints, and all of the priority 2 outlets are tested to see which one
is best used for meeting the target release temperature.

In this example, the priority 1 outlet is given a minimum flow constraint, 
specifying that at least 20% of the total release should go through that outlet.
No other outlets have minimum or maximum flow constraints.  However, the lowest
two outlets have maximum head constraints, such that they cannot be used if they
are deeper than 60 meters.  All outlets were given a minimum head constraint, 
such that they cannot be used unless there is at least 2 meters of depth at the
centerline elevation of the outlet.

So, the result is that outlet #6 is blended with an upper outlet during summer
when the goal is to export warmer water from the lake, and then blended with a
lower outlet when the goal is to export cold water late in the season, subject
to the lowest outlets not being more than 60 meters deep.
--------------------------------------------------------------------------------


For more information, contact:
 Stewart Rounds
 U.S. Geological Survey
 Oregon Water Science Center
 2130 SW 5th Avenue
 Portland, OR 97201
 503-251-3280
 sarounds@usgs.gov

or

 Norman Buccola
 U.S. Geological Survey
 Oregon Water Science Center
 2130 SW 5th Avenue
 Portland, OR 97201
 503-251-3245
 nbuccola@usgs.gov
