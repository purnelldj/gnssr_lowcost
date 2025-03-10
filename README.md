# gnssr_lowcost
Analyze low-cost GNSS data (in NMEA 0183 format) for reflectometry using inverse modelling of SNR data.

The codes can be cited as follows:
Purnell, D.J., Gomez, N., Minarik, W. and Porter, D. (2021). Precise water level measurements using low-cost GNSS antenna arrays (submitted to EGU ESurf).
https://esurf.copernicus.org/preprints/esurf-2020-108/

For any questions or suggestions or issues please contact Dave Purnell:
david.purnell@mail.mcgill.ca

To use in Python: see 'example_python.py'

To use in MATLAB: see 'example_matlab.m'

Codes initially written in Matlab and then translated to Python. Python is much quicker and probably easier to use but Matlab codes were used for the paper results. Tried to make codes as similar as possible but there seems to be some fundamental differences between MATLAB and Python in the lomb-scargle periodogram functions as well as the least squares adjustment. The MATLAB one will probably slowly rot and the Python one may get updated occasionally.

Sample data given in the data directory is one day of data (13th September 2020) from the antenna array at Trois-Rivières (see above paper for more information).

Also see 'code_layout.pdf' for some more information about the key codes in matlab and what they do (inputs/outputs).

All codes written by David Purnell except for: 'ecef2lla.m', 'lla2ecef.m' and the codes contained in 'bspline' and 'fresnel'
See the 'readme' files in those directories for more details
