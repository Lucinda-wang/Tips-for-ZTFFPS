# Tips-for-ZTFFPS
It's an open platform for who are interested in transients. Here we provide some tips to know how to submit object and do light curve analysis.

## Quick Start
```python
import numpy as np

sn_path = '/the/way/to/your/ztffps_file.txt'
MJD  = SN_reader_ztf(sn_path)['mjd']
flux = SN_reader_ztf(sn_path)['flux']

```
