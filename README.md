# gnssr_lowcost
Analyze low-cost GNSS data (in NMEA 0183 format) for reflectometry using inverse modelling.

The codes can be cited as follows:
Purnell, D.J., Gomez, N., Minarik, W. and Porter, D. (2021). Precise water level measurements using low-cost GNSS antenna arrays (submitted to EGU ESurf).

For any questions or suggestions or issues please contact Dave Purnell:
david.purnell@mail.mcgill.ca

See the code 'example_analysis.m' for a step-by-step guide on converting low-cost GNSS data to water level measurements and comparing with measurements from a co-located tide gauge

Also see 'code_layout.pdf' for some more information about the key codes and what they do (inputs/outputs)

All codes written by David Purnell except for: 'ecef2lla.m', 'lla2ecef.m' and the codes contained in 'bspline' and 'fresnel'
See the 'readme' files in those directories for more details
