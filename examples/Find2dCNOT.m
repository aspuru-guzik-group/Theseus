(*
    Conceptual understanding through efficient inverse-design of quantum optical experiments
    
    Example: two-qubit CNOT
    Mario Krenn (mario.krenn@utoronto.ca), Jakob Kottmann, Nora Tischler, Alán Aspuru-Guzik (alan@aspuru.com)
    https://arxiv.org/abs/2005.06443
    mariokrenn.wordpress.com/ https://www.matter.toronto.edu/

    We show an example of Theseus, a highly-efficient inverse-design algorithm for
    quantum optical experiments. In this example, we let Theseus find experimental
    setups for a 2-qubit CNOT gate, using two ancilla gates (related to the experiment
    by Gasparoni et al., Phys. Rev. Lett. 93, 020504). The setups is topologically
    optimized, and uses virtual vertices. The result is a graph which can be translated
    into several different schemes of quantum optics. Here, we dont use topological
    optimization -- see other examples for that. For more infos, see paper.

    Theseus will write a report into a log-file and in the mathematica notebook.
*)

PrintLog[str_]:=(
  Print[str];
  (* log-file *)
  CurrentFileName=NotebookDirectory[]<>"FindHeralded3dimBell_Log.txt"; 
  hFile=OpenAppend[CurrentFileName];
  WriteString[hFile,str<>"\n"];
  Close[hFile];
);

NumberOfSinglePhoton=2;
NumberOfAdditionalPaths=4;
IndividualPaths=Range[NumberOfAdditionalPaths]/.{1->a,2->b,3->c,4->d};
AppendTo[IndividualPaths,Va];
AppendTo[IndividualPaths,Vb];

PossibleComb=Subsets[Range[Length[IndividualPaths]],{2}];
PosCombABCD[n_]:=PossibleComb[[n]]/.{1->a,2->b,3->c,4->d,5->Va,6->Vb,7->Vc,8->Vd}


(*Create Graphs*)

(*First, create edges between a,b,c,d -- afterwards add edges to virtual vertices*)

AllPaths={};
AllPathsWeighted={};
wcount=1;
For[ii=1,ii<NumberOfAdditionalPaths,ii++,
  For[jj=ii+1,jj<=NumberOfAdditionalPaths,jj++,
    CurrPath=ww[wcount++]*IndividualPaths[[ii]][0]*IndividualPaths[[jj]][0];
    AppendTo[AllPathsWeighted,CurrPath];
  ];
];


For[jj=NumberOfSinglePhoton+1,jj<=NumberOfAdditionalPaths,jj++,
  CurrPath=ww[wcount++]a[1]*IndividualPaths[[jj]][0];
  AppendTo[AllPathsWeighted,CurrPath];
  
  CurrPath=ww[wcount++]b[1]*IndividualPaths[[jj]][0];
  AppendTo[AllPathsWeighted,CurrPath];
];

AppendTo[AllPathsWeighted,ww[wcount++]a[0]*b[1]];
AppendTo[AllPathsWeighted,ww[wcount++]a[1]*b[0]];
AppendTo[AllPathsWeighted,ww[wcount++]a[1]*b[1]];


(*Adding edges from vitrual vertices. Those are restricted by SU(2) transformations*)

AllPathsWeighteda0={};
AllPathsWeighteda1={};

AllPathsWeightedb0={};
AllPathsWeightedb1={};
AllPathsWeightedb2={};
For[ii=1,ii<=NumberOfSinglePhoton,ii++,
  CurrPath=ww[wcount]Va[0]*IndividualPaths[[ii]][0];
  AppendTo[AllPathsWeighteda0,CurrPath];
  CurrPath=ww[wcount+1]Va[0]*IndividualPaths[[ii]][1];
  AppendTo[AllPathsWeighteda0,CurrPath];
  
  CurrPath=-ww[wcount+1]*Va[1]*IndividualPaths[[ii]][0];
  AppendTo[AllPathsWeighteda1,CurrPath];
  CurrPath=ww[wcount]*Va[1]*IndividualPaths[[ii]][1];
  AppendTo[AllPathsWeighteda1,CurrPath];
  
  wcount+=2;
  
  CurrPath=ww[wcount]Vb[0]*IndividualPaths[[ii]][0];
  AppendTo[AllPathsWeighteda0,CurrPath];
  CurrPath=ww[wcount+1]Vb[0]*IndividualPaths[[ii]][1];
  AppendTo[AllPathsWeighteda0,CurrPath];
  
  CurrPath=-ww[wcount+1]*Vb[1]*IndividualPaths[[ii]][0];
  AppendTo[AllPathsWeighteda1,CurrPath];
  CurrPath=ww[wcount]*Vb[1]*IndividualPaths[[ii]][1];
  AppendTo[AllPathsWeighteda1,CurrPath];
  
  wcount+=2;
];

For[ii=NumberOfSinglePhoton+1,ii<=NumberOfAdditionalPaths,ii++,
  CurrPath=ww[wcount]Va[0]*IndividualPaths[[ii]][0];
  AppendTo[AllPathsWeighteda0,CurrPath];
  CurrPath=ww[wcount+1]Va[1]*IndividualPaths[[ii]][0];
  AppendTo[AllPathsWeighteda1,CurrPath];
  
  CurrPath=ww[wcount+2]Vb[0]*IndividualPaths[[ii]][0];
  AppendTo[AllPathsWeightedb0,CurrPath];
  CurrPath=ww[wcount+3]Vb[1]*IndividualPaths[[ii]][0];
  AppendTo[AllPathsWeightedb1,CurrPath];
  wcount+=4;
];
wcount-=1;

(*Define the four complete graphs corresponding to 00,01,10,11*)
AllPathsWeighted00=Flatten[{AllPathsWeighted,AllPathsWeighteda0,AllPathsWeightedb0}];
AllPathsWeighted01=Flatten[{AllPathsWeighted,AllPathsWeighteda0,AllPathsWeightedb1}];
AllPathsWeighted10=Flatten[{AllPathsWeighted,AllPathsWeighteda1,AllPathsWeightedb0}];
AllPathsWeighted11=Flatten[{AllPathsWeighted,AllPathsWeighteda1,AllPathsWeightedb1}];


TriggerNoMultiPhoton={Va[_]^n_->0,Vb[_]^n_->0,c[_]^n_->0,d[_]^n_->0,c[n1_]*c[n2_]->0,d[n1_]*d[n2_]->0};

(*generate the higher-order terms of graph weight function \[CapitalPhi](\[Omega]) *)
(*To speed up, we already remove terms that cannot make contributions, such as double excitation of the virtual vertices (they correspond to exactly one incoming photon*)

FS00=Expand[Total[AllPathsWeighted00]];
FS00e2=Expand[FS00*FS00]/.TriggerNoMultiPhoton;
FS00e3=Expand[FS00e2*FS00]/.TriggerNoMultiPhoton;
FS00e4=Expand[FS00e3*v[0]]/.TriggerNoMultiPhoton;
TriggerableState00=Expand[FS00e4*ZERO]/.{ZERO*c[0]*d[0]*Va[dVa_]*Vb[dVb_]v[0]->c[0]*d[0]*Va[dVa]*Vb[dVb]*v[0]}/.{ZERO->0};
Tmp1=TriggerableState00/.{Plus->List}/.{ww[_]->1};
NormTmp=Tmp1/.{a[_]->1,b[_]->1,c[_]->1,d[_]->1,e[_]->1,f[_]->1,Va[_]->1,Vb[_]->1,v[0]->1};
AllABCCombinations00=DeleteDuplicates[Tmp1/NormTmp];
Print["Finished 00"];

FS01=Expand[Total[AllPathsWeighted01]];
FS01e2=Expand[FS01*FS01]/.TriggerNoMultiPhoton;
FS01e3=Expand[FS01e2*FS01]/.TriggerNoMultiPhoton;
FS01e4=Expand[FS01e3*v[0]]/.TriggerNoMultiPhoton;
TriggerableState01=Expand[FS01e4*ZERO]/.{ZERO*c[0]*d[0]*Va[dVa_]*Vb[dVb_]v[0]->c[0]*d[0]*Va[dVa]*Vb[dVb]*v[0]}/.{ZERO->0};
Tmp1=TriggerableState01/.{Plus->List}/.{ww[_]->1};
NormTmp=Tmp1/.{a[_]->1,b[_]->1,c[_]->1,d[_]->1,e[_]->1,f[_]->1,Va[_]->1,Vb[_]->1,v[0]->1};
AllABCCombinations01=DeleteDuplicates[Tmp1/NormTmp];
Print["Finished 01"];

FS10=Expand[Total[AllPathsWeighted10]];
FS10e2=Expand[FS10*FS10]/.TriggerNoMultiPhoton;
FS10e3=Expand[FS10e2*FS10]/.TriggerNoMultiPhoton;
FS10e4=Expand[FS10e3*v[0]]/.TriggerNoMultiPhoton;
TriggerableState10=Expand[FS10e4*ZERO]/.{ZERO*c[0]*d[0]*Va[dVa_]*Vb[dVb_]v[0]->c[0]*d[0]*Va[dVa]*Vb[dVb]*v[0]}/.{ZERO->0};
Tmp1=TriggerableState10/.{Plus->List}/.{ww[_]->1};
NormTmp=Tmp1/.{a[_]->1,b[_]->1,c[_]->1,d[_]->1,e[_]->1,f[_]->1,Va[_]->1,Vb[_]->1,v[0]->1};
AllABCCombinations10=DeleteDuplicates[Tmp1/NormTmp];
Print["Finished 10"];

FS11=Expand[Total[AllPathsWeighted11]];
FS11e2=Expand[FS11*FS11]/.TriggerNoMultiPhoton;
FS11e3=Expand[FS11e2*FS11]/.TriggerNoMultiPhoton;
FS11e4=Expand[FS11e3*v[0]]/.TriggerNoMultiPhoton;
TriggerableState11=Expand[FS11e4*ZERO]/.{ZERO*c[0]*d[0]*Va[dVa_]*Vb[dVb_]v[0]->c[0]*d[0]*Va[dVa]*Vb[dVb]*v[0]}/.{ZERO->0};
Tmp1=TriggerableState11/.{Plus->List}/.{ww[_]->1};
NormTmp=Tmp1/.{a[_]->1,b[_]->1,c[_]->1,d[_]->1,e[_]->1,f[_]->1,Va[_]->1,Vb[_]->1,v[0]->1};
AllABCCombinations11=DeleteDuplicates[Tmp1/NormTmp];
Print["Finished 11"];


(*For each graph, create all available terms from the graph weight function, including the target term. those will be used for the fidelity later*)

AllEquations={};
TargetEquations={};
ExpandedStateZERO00=Expand[TriggerableState00*ZERO];
AllABCCombinations=AllABCCombinations00;
For[ii=1,ii<=Length[AllABCCombinations],ii++,
  newEq=ExpandedStateZERO00/.{ZERO*AllABCCombinations[[ii]]->1}/.{ZERO->0};
  If[(AllABCCombinations[[ii]]==a[0]*b[0]*c[0]*d[0]*v[0]*Va[0]*Vb[0]),
    AppendTo[TargetEquations,newEq];
    Print["x00"]
  ];
  AppendTo[AllEquations,newEq];
];
AllEquations00=AllEquations/.{a[_]->0,b[_]->0,c[_]->0,d[_]->0,e[_]->0,f[_]->0,g[_]->0,h[_]->0,i[_]->0,j[_]->0,a[_]^n_->0,b[_]^n_->0,c[_]^n_->0,d[_]^n_->0,e[_]^n_->0,f[_]^n_->0,g[_]^n_->0,h[_]^n_->0,i[_]^n_->0,j[_]^n_->0};
TargetEquations00=TargetEquations/.{a[_]->0,b[_]->0,c[_]->0,d[_]->0,e[_]->0,f[_]->0,g[_]->0,h[_]->0,i[_]->0,j[_]->0,a[_]^n_->0,b[_]^n_->0,c[_]^n_->0,d[_]^n_->0,e[_]^n_->0,f[_]^n_->0,g[_]^n_->0,h[_]^n_->0,i[_]^n_->0,j[_]^n_->0};




AllEquations={};
TargetEquations={};
ExpandedStateZERO01=Expand[TriggerableState01*ZERO];
AllABCCombinations=AllABCCombinations01;
For[ii=1,ii<=Length[AllABCCombinations],ii++,
  newEq=ExpandedStateZERO01/.{ZERO*AllABCCombinations[[ii]]->1}/.{ZERO->0};
  If[(AllABCCombinations[[ii]]==a[0]*b[1]*c[0]*d[0]*v[0]*Va[0]*Vb[1]),
    AppendTo[TargetEquations,newEq];
    Print["x01"];
  ];
  AppendTo[AllEquations,newEq];
];
AllEquations01=AllEquations/.{a[_]->0,b[_]->0,c[_]->0,d[_]->0,e[_]->0,f[_]->0,g[_]->0,h[_]->0,i[_]->0,j[_]->0,a[_]^n_->0,b[_]^n_->0,c[_]^n_->0,d[_]^n_->0,e[_]^n_->0,f[_]^n_->0,g[_]^n_->0,h[_]^n_->0,i[_]^n_->0,j[_]^n_->0};
TargetEquations01=TargetEquations/.{a[_]->0,b[_]->0,c[_]->0,d[_]->0,e[_]->0,f[_]->0,g[_]->0,h[_]->0,i[_]->0,j[_]->0,a[_]^n_->0,b[_]^n_->0,c[_]^n_->0,d[_]^n_->0,e[_]^n_->0,f[_]^n_->0,g[_]^n_->0,h[_]^n_->0,i[_]^n_->0,j[_]^n_->0};



AllEquations={};
TargetEquations={};
ExpandedStateZERO10=Expand[TriggerableState10*ZERO];
AllABCCombinations=AllABCCombinations10;
For[ii=1,ii<=Length[AllABCCombinations],ii++,
  newEq=ExpandedStateZERO10/.{ZERO*AllABCCombinations[[ii]]->1}/.{ZERO->0};
  If[(AllABCCombinations[[ii]]==a[1]*b[1]*c[0]*d[0]*v[0]*Va[1]*Vb[0]),
    AppendTo[TargetEquations,newEq];
    Print["x10"];
  ];
  AppendTo[AllEquations,newEq];
];
AllEquations10=AllEquations/.{a[_]->0,b[_]->0,c[_]->0,d[_]->0,e[_]->0,f[_]->0,g[_]->0,h[_]->0,i[_]->0,j[_]->0,a[_]^n_->0,b[_]^n_->0,c[_]^n_->0,d[_]^n_->0,e[_]^n_->0,f[_]^n_->0,g[_]^n_->0,h[_]^n_->0,i[_]^n_->0,j[_]^n_->0};
TargetEquations10=TargetEquations/.{a[_]->0,b[_]->0,c[_]->0,d[_]->0,e[_]->0,f[_]->0,g[_]->0,h[_]->0,i[_]->0,j[_]->0,a[_]^n_->0,b[_]^n_->0,c[_]^n_->0,d[_]^n_->0,e[_]^n_->0,f[_]^n_->0,g[_]^n_->0,h[_]^n_->0,i[_]^n_->0,j[_]^n_->0};



AllEquations={};
TargetEquations={};
ExpandedStateZERO11=Expand[TriggerableState11*ZERO];
AllABCCombinations=AllABCCombinations11;
For[ii=1,ii<=Length[AllABCCombinations],ii++,
  newEq=ExpandedStateZERO11/.{ZERO*AllABCCombinations[[ii]]->1}/.{ZERO->0};
  If[(AllABCCombinations[[ii]]==a[1]*b[0]*c[0]*d[0]*v[0]*Va[1]*Vb[1]),
    AppendTo[TargetEquations,newEq];
    Print["x11"];
  ];
  AppendTo[AllEquations,newEq];
];
AllEquations11=AllEquations/.{a[_]->0,b[_]->0,c[_]->0,d[_]->0,e[_]->0,f[_]->0,g[_]->0,h[_]->0,i[_]->0,j[_]->0,a[_]^n_->0,b[_]^n_->0,c[_]^n_->0,d[_]^n_->0,e[_]^n_->0,f[_]^n_->0,g[_]^n_->0,h[_]^n_->0,i[_]^n_->0,j[_]^n_->0};
TargetEquations11=TargetEquations/.{a[_]->0,b[_]->0,c[_]->0,d[_]->0,e[_]->0,f[_]->0,g[_]->0,h[_]->0,i[_]->0,j[_]->0,a[_]^n_->0,b[_]^n_->0,c[_]^n_->0,d[_]^n_->0,e[_]^n_->0,f[_]^n_->0,g[_]^n_->0,h[_]^n_->0,i[_]^n_->0,j[_]^n_->0};


Prefactor=1/4*(1/Total[Abs[Flatten[{AllEquations00,AllEquations01,AllEquations10,AllEquations11}]]^2]);
Overlap=Total[Flatten[{TargetEquations00,TargetEquations01,TargetEquations10,TargetEquations11}]];


(*The fidelity of the whole process is written as combination of the fidelity of the individual graphs.*)

TotalFidelity=Prefactor*Abs[Overlap]^2;
Loss=(1-TotalFidelity);

vars=Sort[DeleteDuplicates[Cases[Loss,_ww,{0,Infinity}]]];
AllStandardVars=vars[[1;;13]];

TotalSumVars=Total[Abs[AllStandardVars]];

alpha=0.005;
Loss2=(1-TotalFidelity)+alpha*TotalSumVars;

LenOfVariables=Length[vars];

Off[FindMinimum::lstol];
Off[FindMinimum::nrnum];
Off[FindMinimum::fmgz];
Off[FindMinimum::cvmit];
StartTime=AbsoluteTime[];

Print["Entering Main Loop"];
res=0;
cursum=0;
NumOfRepeatInit=250;(*Maximally 250 iterations when attempting to remove edge*)

ZeroList={};
NumOfRepeat=NumOfRepeatInit;
While[NumOfRepeat>0,
  NumOfRepeat-=1;
  RemoveVar=RandomChoice[vars];
  
  (*We chose a random edge that will be removed. could be chosen via the results of previous optimization for speedup. See heralded-3dBell example. *)
  CurrZeroList=Append[ZeroList,RemoveVar->0];
  newLoss2=Loss2/.CurrZeroList;
  newLoss=Loss/.CurrZeroList;
  
  newvars=Sort[DeleteDuplicates[Cases[newLoss2,_ww,{0,Infinity}]]];
  
  SmallVars=0;
  For[ii=1,ii<=Length[newvars],ii++,
    If[newvars[[ii,1]]<=13,SmallVars++;];
  ];
  BigVars=Length[newvars]-SmallVars;
  rndVars=Flatten[{RandomReal[{-0.1,0.1},SmallVars],RandomReal[{-1,1},Length[newvars]-SmallVars]}];
  fullvars=Transpose[{newvars,rndVars}];
  
  newsol=FindMinimum[newLoss2,fullvars,AccuracyGoal->4,PrecisionGoal->4];
  newres=newLoss/.newsol[[2]];
  cursum=TotalSumVars/.CurrZeroList/.newsol[[2]];
  If[Mod[NumOfRepeat,25]==0,
    Print["|\[Omega]Subscript[|, 1]: "<>ToString[cursum,InputForm]<>"; Fidelity="<>ToString[1-newres,InputForm]];
  ];
  If[cursum<1&&newres<0.05,
    ZeroList=CurrZeroList;
    newgoodloss=newLoss;
    newgoodloss2=newLoss2;
    newgoodsol=newres;
    vars=newvars;
    Print[ToString[Length[vars]]<>" edges ("<>ToString[NumOfRepeat]<>"), edges from Va/Vb ("<>ToString[BigVars]<>"): "<>ToString[cursum,InputForm]<>"; Fidelity="<>ToString[1-newres,InputForm]];
    NumOfRepeat=NumOfRepeatInit;
    BigVars=0;
    For[ii=1,ii<=Length[newvars],ii++,
      If[newvars[[ii,1]]>13,BigVars++;];
    ];
    CurrentFileName=NotebookDirectory[]<>"Find2dCNOT_Log.txt"; 
    hFile=OpenAppend[CurrentFileName];
    WriteString[hFile,"#(edges)=: "<>ToString[Length[vars]]<>"\n"];
    WriteString[hFile,"#(virtualEdges)=: "<>ToString[BigVars]<>"\n"];
    WriteString[hFile,"ZeroList: "<>ToString[ZeroList]<>"\n"];
    WriteString[hFile,"newgoodsol: "<>ToString[newsol[[2]],InputForm]<>"\n\n\n\n\n"];
    Close[hFile];
  
  ];
];

EndTime=AbsoluteTime[];
Print["Time: "<>ToString[EndTime-StartTime]<>" sec."];
Print["Final: "<>ToString[LenOfVariables-Length[ZeroList]]<>"/"<>ToString[LenOfVariables]<>" variables."];
