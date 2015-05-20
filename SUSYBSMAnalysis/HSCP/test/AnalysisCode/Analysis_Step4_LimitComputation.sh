root -l -b << EOF
   TString makeshared(gSystem->GetMakeSharedLib());
   TString dummy = makeshared.ReplaceAll("-W ", "-Wno-deprecated-declarations -Wno-deprecated -D__USE_XOPEN2K8 ");
   TString dummy = makeshared.ReplaceAll("-Wshadow ", " -std=c++0x ");
   cout << "Compilling with the following arguments: " << makeshared << endl;
   gSystem->SetMakeSharedLib(makeshared);
   gSystem->SetIncludePath( "-I$ROOFITSYS/include" );
  .x Analysis_Step4_LimitComputation.C++("Final7TeV", "", "");
  .x Analysis_Step4_LimitComputation.C++("Final8TeV", "", "");
  .x Analysis_Step4_LimitComputation.C++("FinalCOMB", "", "");
EOF

