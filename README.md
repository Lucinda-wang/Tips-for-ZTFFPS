# Tips-for-ZTFFPS
It's an open platform for who are interested in transients. Here we provide some tips to know how to submit object and do light curve analysis.

## Description
Before you start your program, you will first need to register by email and indicate your request for registration. (Remark: Strongly advised you register with your institutional email address)
1. Register - Access to ztf@ipac.caltech.edu
2. Receive - After your registration got accept, you will also got your own **password** and **submission account**
3. Double check - Ensure your target can be observed by ZTF based on the coordinate.
4. Record - Record the following parameters before submission.
   - RA: Right Ascension of the target position. (decimal degrees; J2000)
   - DEC: Declination of the target position. (decimal degrees; J2000)
   - JD_start: Starting date for observation in unit of Julian Date. (decimal days)
   - JD_end Ending date for observation in unit of Julian Date. (decimal days)
   - Email_address: The email address that you registered with.
   - Userpass: The password you received

## It's time to submit HOHOHO~
Here is the graphical user interface (GUI) form for you to submit your request:
https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi

After you enter the service, you need to write the first specific username and password.
username: ztffps

password: dontgocrazy!

## Check your HOHOHO~
If you would like to check if your target had been submitted already, check this out:
https://ztfweb.ipac.caltech.edu/cgi-bin/getForcedPhotometryRequests.cgi


## Quick Start
```python
import numpy as np

sn_path = '/the/way/to/your/ztffps_file.txt'
MJD  = SN_reader_ztf(sn_path)['mjd']
flux = SN_reader_ztf(sn_path)['flux']

```

## Links
- Reference: https://ui.adsabs.harvard.edu/abs/2019PASP..131a8003M/abstract
