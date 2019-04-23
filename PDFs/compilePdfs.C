void compilePdfs()
{
	cout<<"===========	First Line 		================"<<endl;
	gSystem->AddIncludePath("-I$ROOFITSYS/include");
	cout<<"\n\n===========	processing HWWLVJRooPdfs.cxx	========"<<endl;
	gROOT->ProcessLine(".L HWWLVJRooPdfs.cxx+");
	cout<<"\n\n===========	processing PdfDiagonalizer.cc	========"<<endl;
	gROOT->ProcessLine(".L PdfDiagonalizer.cc+");
	cout<<"\n\n===========	processing Util.cxx		========"<<endl;
	gROOT->ProcessLine(".L Util.cxx+");
}
