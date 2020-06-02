void example()
{
  for(int i=0;i<nFiles;i++) { // 1 file = 1 bunch
    treeMain=file[i]->GetObject("tree");

    //// take tail component from other files
    for (int j=0;j<nFiles-1;j++) {
      k=i+1;
      if(k>nFiles) k-=nFiles;
      treeTail=file[k]->GetObject("tree");
      ///
      tLimit[2] = {1170*(j+1),1170*(j+2)};
    }
  }
}
