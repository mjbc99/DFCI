### Description of scripts

# sup.moa2.R
# Author: Azfar
# Description: ????
# Sourced by Makena_Enrichment.R,Makena_Enrichment_no-normal.R


# global-test.R
# Author:Makena
# Description: Run GSA correcting for %Tumor Purity
# Depreciated Use. global-test-2.R


# global-test-2.R
# Author:Makena
# Description: Run GSA correcting for %Tumor Purity.  Except changed basal to lum#

# Makena_Enrichment.R
# Author: Makena and Azfar
# Description: Run GSA using MsigDB/ImmunoDB, GeneSigDB and BindeaDB on LumA IDC v ILC data. Includes quantiile normalization step
# Depends on: sup.moa2.R

# Makena_Enrichment_no-normal.R
# Author: Makena and Azfar
# Description: Run GSA using MsigDB/ImmunoDB, GeneSigDB and BindeaDB on LumA IDC v ILC data.  Same as Makena_Enrichment.R but without quantiile normalization step
# Depends on: sup.moa2.R