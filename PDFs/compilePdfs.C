void compilePdfs()
{
	cout<<"===========	First Line 		================"<<endl;
	gSystem->AddIncludePath("-I$ROOFITSYS/include");
	cout<<"===========	processing HWWLVJRooPdfs.cxx	========"<<endl;
	gROOT->ProcessLine(".L HWWLVJRooPdfs.cxx+");
	cout<<"===========	processing PdfDiagonalizer.cc	========"<<endl;
	gROOT->ProcessLine(".L PdfDiagonalizer.cc+");
	cout<<"===========	processing Util.cxx		========"<<endl;
	gROOT->ProcessLine(".L Util.cxx+");
}
