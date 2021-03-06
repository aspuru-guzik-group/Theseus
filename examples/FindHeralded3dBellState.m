(*
    Conceptual understanding through efficient inverse-design of quantum optical experiments
    
    Example: Heralded 3-dimensional Bell state
    Mario Krenn (mario.krenn@utoronto.ca), Jakob Kottmann, Nora Tischler, Alán Aspuru-Guzik (alan@aspuru.com)
    https://arxiv.org/abs/2005.06443
    mariokrenn.wordpress.com/ https://www.matter.toronto.edu/
    We show an example of Theseus, a highly-efficient inverse-design algorithm for quantum optical experiments.
    In thisexample, we let Theseus find experimental setups a heralded, 3-dimensional Bell state in path a,b --
    using six ancilla photons. The result is a graph which can be translated into several different schemes of
    quantum optics. Here, we dont use topological optimization -- see other examples for that.
    For more infos, see paper.
    
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

PossibleComb=Subsets[Range[8],{2}];
PosCombABCD[n_]:=PossibleComb[[n]]/.{1->a,2->b,3->c,4->d,5->e,6->f,7->g,8->h}

(*build Graph *)
AllPaths={};
For[ii=1,ii<=Length[PossibleComb],ii++,
  For[d1=0,d1<=2,d1++,
    For[d2=0,d2<=2,d2++,
      AppendTo[AllPaths,PosCombABCD[ii][[1]][d1]*PosCombABCD[ii][[2]][d2]];
    ];
  ];
];

(*vertices c,d,e,f,g,h only need one mode number, as they are the heralding detectors*)
AllPaths=DeleteDuplicates[AllPaths/.{c[_]->c[0],d[_]->d[0],e[_]->e[0],f[_]->f[0],g[_]->g[0],h[_]->h[0]}];

(*give weights to edges*)
AllPathsWeighted={};
VarList={};
For[ii=1,ii<=Length[AllPaths],ii++,
  AppendTo[AllPathsWeighted,ww[ii]*AllPaths[[ii]]];
  AppendTo[VarList,ww[ii]];
];


(* Create weight function phi(Omega), up to 4th order *)
AllEdgesSum=Total[AllPathsWeighted]; Print["Created 1st order"];
DoublePaths=Expand[AllEdgesSum^2]; Print["Created 2nd order"];
TriplePaths=Expand[AllEdgesSum^3]; Print["Created 3rd order"];
QuadPaths=Expand[AllEdgesSum^4]; Print["Created 4th order"];
FullState=Expand[(AllEdgesSum+DoublePaths/2+TriplePaths/6+QuadPaths/24)*v[0]]; Print["Created full weight function phi(Omega)"]


(*Condition on a click in the heralding detectors*)
TriggerableStateC=Expand[FullState*ZEROc]/.{ZEROc*c[0]*v[0]->v[0],ZEROc*c[0]^n_*v[0]->c[0]^(n-1)*v[0]}/.{ZEROc->0}; PrintLog["Detector C heralded"]
TriggerableStateD=Expand[TriggerableStateC*ZEROd]/.{ZEROd*d[0]*v[0]->v[0],ZEROd*d[0]^n_*v[0]->d[0]^(n-1)*v[0]}/.{ZEROd->0}; PrintLog["Detector D heralded"]
TriggerableStateE=Expand[TriggerableStateD*ZEROe]/.{ZEROe*e[0]*v[0]->v[0],ZEROe*e[0]^n_*v[0]->e[0]^(n-1)*v[0]}/.{ZEROe->0}; PrintLog["Detector E heralded"]
TriggerableStateF=Expand[TriggerableStateE*ZEROf]/.{ZEROf*f[0]*v[0]->v[0],ZEROf*f[0]^n_*v[0]->f[0]^(n-1)*v[0]}/.{ZEROf->0}; PrintLog["Detector F heralded"]
TriggerableStateG=Expand[TriggerableStateF*ZEROg]/.{ZEROg*g[0]*v[0]->v[0],ZEROg*g[0]^n_*v[0]->g[0]^(n-1)*v[0]}/.{ZEROg->0}; PrintLog["Detector G heralded"]
TriggerableStateH=Expand[TriggerableStateG*ZEROh]/.{ZEROh*h[0]*v[0]->v[0],ZEROh*h[0]^n_*v[0]->h[0]^(n-1)*v[0]}/.{ZEROh->0}; PrintLog["Detector H heralded"]
TriggerableState=TriggerableStateH;


(*Find all possible terms in quantum state*)
AllABCCombinations=DeleteDuplicates[TriggerableState/.{ww[_]->1}/.{Plus->List}];
NormCombTmp=AllABCCombinations/.{a[_]->1,b[_]->1,c[_]->1,d[_]->1,e[_]->1,f[_]->1,g[_]->1,h[_]->1,v[_]->1};
NormComb=ConstantArray[1,Length[AllABCCombinations]]/NormCombTmp;
AllABCCombinations=AllABCCombinations*NormComb;
PrintLog["Number of possible terms: "<>ToString[Length[AllABCCombinations]]];
PrintLog["Number of edges: "<>ToString[Length[VarList]]];


AllEquations={};
TargetEquations={};
ExpandedStateZERO=Expand[TriggerableState*ZERO];

(*Create Lists of terms that are used in fidelity*)
(*in particular, find terms of 3d Bell state in terms of graph weights*)
For[ii=1, ii<=Length[AllABCCombinations], ii++,
  newEq=ExpandedStateZERO/.{ZERO*AllABCCombinations[[ii]]->1}/.{ZERO->0};
  If[(AllABCCombinations[[ii]]==a[0]*b[0]*v[0])||(AllABCCombinations[[ii]]==a[1]*b[1]*v[0])||(AllABCCombinations[[ii]]==a[2]*b[2]*v[0]),
    AppendTo[TargetEquations,newEq];
  ];
  AppendTo[AllEquations,newEq];
];
AllEquations=AllEquations/.{a[_]->0,b[_]->0,c[_]->0,d[_]->0,e[_]->0,f[_]->0,g[_]->0,h[_]->0,i[_]->0,j[_]->0,a[_]^n_->0,b[_]^n_->0,c[_]^n_->0,d[_]^n_->0,e[_]^n_->0,f[_]^n_->0,g[_]^n_->0,h[_]^n_->0,i[_]^n_->0,j[_]^n_->0};
TargetEquations=TargetEquations/.{a[_]->0,b[_]->0,c[_]->0,d[_]->0,e[_]->0,f[_]->0,g[_]->0,h[_]->0,i[_]->0,j[_]->0,a[_]^n_->0,b[_]^n_->0,c[_]^n_->0,d[_]^n_->0,e[_]^n_->0,f[_]^n_->0,g[_]^n_->0,h[_]^n_->0,i[_]^n_->0,j[_]^n_->0};




NormalisationConstant=Total[Abs[AllEquations]^2];
Fidelity=Total[Abs[TargetEquations]]^2/(Length[TargetEquations]*NormalisationConstant);
Loss=(1-Fidelity);
vars=Sort[DeleteDuplicates[Cases[NormalisationConstant,_ww,{0,Infinity}]]];
TotalSumVars=Total[Abs[vars]];

alpha=0.5;
Loss2=(1-Fidelity)+alpha*TotalSumVars; (* Loss contains L1 regularisation *)


LenOfVariables=Length[vars];

Off[FindMinimum::lstol];
Off[FindMinimum::nrnum];
Off[FindMinimum::fmgz];
Off[FindMinimum::cvmit];
StartTime=AbsoluteTime[];

Print["Entering Main Loop"];
ZeroList={}; (*Collection of edges that are removed*)
NumOfRepeatInitial=500;

OptimizedEdgeWeights=ConstantArray[0,Length[vars]];
NumOfRepeat=NumOfRepeatInitial;
While[NumOfRepeat>0,
  idx=(NumOfRepeatInitial-NumOfRepeat)+1;
  If[((NumOfRepeatInitial-NumOfRepeat)<Floor[Length[vars]/2])&&RandomReal[]<0.75,
    (*Here we remove a n edge according to results from previous optimization*)
    pos=Ordering[OptimizedEdgeWeights];
    RemoveVar=vars[[pos[[idx]]]];
    ,
    (*Here we remove a random edge.*)
    RemoveVar=RandomChoice[vars];
  ];

  NumOfRepeat-=1;

  CurrZeroList=Append[ZeroList,RemoveVar->0];
  newLoss2=Loss2/.CurrZeroList;
  newLoss=Loss/.CurrZeroList;

  newvars=Sort[DeleteDuplicates[Cases[newLoss2,_ww,{0,Infinity}]]];
  fullvars=Transpose[{newvars,RandomReal[{-0.05,0.05},Length[newvars]]}];
  newsol=FindMinimum[newLoss2,fullvars,AccuracyGoal->4,PrecisionGoal->4];(*quasi-newton iterative method*)

  newres=newLoss/.newsol[[2]];
  cursum=TotalSumVars/.CurrZeroList/.newsol[[2]];
  If[Mod[NumOfRepeat,25]==0,
    PrintLog[ToString[Length[vars]]<>" edges ("<>ToString[NumOfRepeat]<>" more trials): |\[Omega]Subscript[|, 1]="<>ToString[cursum,InputForm]<>"; Fidelity="<>ToString[1-newres,InputForm]];
  ];

  If[cursum<1&&newres<0.05,
    ZeroList=CurrZeroList;
    vars=newvars;
    OptimizedEdgeWeights=Abs[vars/.newsol[[2]]];
    PrintLog["Removed one edge, new #(edges)="<>ToString[Length[vars]]<>"; |\[Omega]Subscript[|, 1]="<>ToString[cursum,InputForm]<>"; Fidelity="<>ToString[1-newres,InputForm]];
    NumOfRepeat=NumOfRepeatInitial;

    CurrentFileName=NotebookDirectory[]<>"FindHeralded3dimBell_Log.txt"; 
    hFile=OpenAppend[CurrentFileName];
    WriteString[hFile,"#(edges)=: "<>ToString[Length[vars]]<>"\n"];
    WriteString[hFile,"ZeroList: "<>ToString[ZeroList]<>"\n"];
    WriteString[hFile,"newgoodsol: "<>ToString[newsol[[2]],InputForm]<>"\n\n\n\n\n"];
    Close[hFile];
  ];
];

PrintLog["Time: "<>ToString[AbsoluteTime[]-StartTime]<>" sec."];
PrintLog["Final: "<>ToString[LenOfVariables-Length[ZeroList]]<>"/"<>ToString[LenOfVariables]<>" variables."];
