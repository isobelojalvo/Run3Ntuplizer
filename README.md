Installation:

```

cmsrel CMSSW_11_2_0
cd CMSSW_11_2_0/src/
cmsenv
git cms-init
git remote add pallabidas https://github.com/pallabidas/cmssw.git
git cms-merge-topic pallabidas:l1t-integration-test-17May
cd L1Trigger
git clone -b boostedTau git@github.com:isobelojalvo/Run3Ntuplizer.git
cd ..
scram b -j 12

cd L1Trigger/Run3Ntuplizer/test

## to run zero bias data
cmsRun testL1TCaloSummary-ZeroBias.py

## to run boosted ggHbb events
cmsRun testL1TCaloSummary-ggHbb.py


```
