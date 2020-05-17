
(*
    Conceptual understanding through efficient inverse-design of quantum optical experiments
    
    Example: Post-selected high-dimensional three-photon entanglement
    Mario Krenn (mario.krenn@utoronto.ca), Jakob Kottmann, Nora Tischler, Alán Aspuru-Guzik (alan@aspuru.com)
    https://arxiv.org/abs/2005.06443
    mariokrenn.wordpress.com/ https://www.matter.toronto.edu/

  We show an example of Theseus, a highly-efficient inverse-design algorithm for
  quantum optical experiments. In this example, we let Theseus find experimental 
  setups for maximally entangled, high-dimensional three-photon quantums states
  (following the definition of Huber&de Vicente, Phys. Rev. Lett. 110, 030501 (2013)).
  The result is a graph which can be translated into several different schemes of
  quantum optics. Here, we dont use topological optimization -- see other examples
  for that. For more infos, see paper.

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



CreateSRVState[localDim_]:=(
  (*
    localDim={1,a,b,c} where a,b,c define the SRV of a state;
    This function returns a maximally entangled quantum state with the given SRV;

    For example:
    CreateSRVState[{1,3,3,2}]={{0,2,0},{1,1,0},{2,0,1}};
    That corresponds to |psi> = |0,2,0> + |1,1,0> + |2,0,1>
  *)
  ccDim=Range[localDim[[3]]]-1;
  ddDim=Range[localDim[[4]]]-1;
  cdtuples=Tuples[{ccDim,ddDim}];
  Bestrs=RandomSample[cdtuples[[1;;localDim[[2]]]]];
  BestDiffSRVnum=2*localDim[[2]];

  rs=Bestrs;
  cc=0;
  IsRunning=True;
  While[IsRunning,
  
  num1=Length[DeleteDuplicates[rs[[;;,1]]]];num2=Length[DeleteDuplicates[rs[[;;,2]]]];
  DiffSRVnum=(localDim[[3]]-num1)+(localDim[[4]]-num2);
  If[DiffSRVnum==0,
    IsRunning=False;
    ,
    If[DiffSRVnum<=BestDiffSRVnum,
      BestDiffSRVnum=DiffSRVnum;
      Bestrs=rs;
    ];

    If[cc>100&RandomReal[]<0.001,
      Bestrs=RandomSample[cdtuples[[1;;localDim[[2]]]]];
      BestDiffSRVnum=2*localDim[[2]];
      Print["CreateSRVState["<>ToString[localDim]<>"]: Restart (cc="<>ToString[cc]<>")"];
      Pause[0.5];
    ];
  
    rs=Bestrs;
    complSets=Complement[cdtuples,rs];
    source=RandomChoice[complSets];
    target=RandomInteger[{1,Length[rs]}];
  
    rs[[target]]=source;
  ];
  
  cc++;
  If[Mod[cc,1000]==0,
    Print["CreateSRVState: cc="<>ToString[cc]];
    Pause[0.5];
  ];
  ];
  AllStatesCSRVS={};
  For[iic=0,iic<localDim[[2]],iic++,
    AppendTo[AllStatesCSRVS,Flatten[{iic,rs[[iic+1]]}]];
  ];
  Return[AllStatesCSRVS];
)

Off[FindMinimum::lstol];
Off[FindMinimum::nrnum];
Off[FindMinimum::fmgz];
Off[FindMinimum::cvmit];

TotalGood={};
TotalTime={};

(* How many iterations are we trying until we give up on state? (Many states are
found within the first iteration) *)
MaxSaveCount=100; 

PrintLog["SRV: 1-Fidelity (time in sec), iterations"];

TotalTimeStart=AbsoluteTime[];
  
NumOfGoodStates=0;
NumOfBadStates=0;
dimBMax=4; (*The largest local dimension we look for is d=4. *)
(*We loop through all possible states with, with restriction d=4 *)

For[dimB=2,dimB<=dimBMax,dimB++,
  For[dimC=2,dimC<=dimB,dimC++,
    For[dimD=2,dimD<=dimC,dimD++,

      StartTime=AbsoluteTime[];
      localDim={1,dimB,dimC,dimD}; (*SRV={dimB,dimC,dimD} with particle A considered as a trigger*)

      (* SRV states with that do not fulfil the following restriction cannot be generated, see Phys.Rev.A 99,032338 (2019) *)
      GraphRestriction1=1+Min[1+(localDim[[2]]-localDim[[3]]),localDim[[4]]]+Min[1+(localDim[[2]]-localDim[[4]]),localDim[[3]]-1]>=localDim[[2]];
      (* SRV states with that do not fulfil the following restriction cannot be written down, see Phys.Rev.Lett.110,030501 (2013) *)
      GraphRestriction2=(localDim[[3]]*localDim[[4]])>=localDim[[2]];

      If[GraphRestriction1&&GraphRestriction2,

        saveCount=0;
        L1val=1; (*Loss=(1-Fidelity) *)
        While[L1val>0.01&&saveCount<=MaxSaveCount,

          AllStates=CreateSRVState[localDim];

          (*Creating the complete graph -- given the restiction of local dimensions *)
          AllPaths={};
          For[da=0,da<localDim[[1]],da++,
            For[db=0,db<localDim[[2]],db++,
              AppendTo[AllPaths,a[da]*b[db]];
            ];
          For[dc=0,dc<localDim[[3]],dc++,
            AppendTo[AllPaths,a[da]*c[dc]];
          ];
          For[dd=0,dd<localDim[[4]],dd++,
            AppendTo[AllPaths,a[da]*d[dd]];
          ];
        ];

        For[db=0,db<localDim[[2]],db++,
          For[dc=0,dc<localDim[[3]],dc++,
            AppendTo[AllPaths,b[db]*c[dc]];
          ];
          For[dd=0,dd<localDim[[4]],dd++,
            AppendTo[AllPaths,b[db]*d[dd]];
          ];
        ];

        For[dc=0,dc<localDim[[3]],dc++,
          For[dd=0,dd<localDim[[4]],dd++,
            AppendTo[AllPaths,c[dc]*d[dd]];
          ];
        ];
        
        (*Every edge gets a weight *)
        AllPathsWeighted={};
        VarList={};
        For[ii=1,ii<=Length[AllPaths],ii++,
          AppendTo[AllPathsWeighted,ww[ii]*AllPaths[[ii]]];
          AppendTo[VarList,ww[ii]];
        ];

        (*The graph weight function Phi(omega) in the paper *)


      
        (*Creating all four-fold detections for this SRV*)
        (*First create the second order state, i.e. two excited edges*)
        DoublePaths=Expand[Total[AllPathsWeighted]^2]/.{a[_]^n_->0,b[_]^n_->0,c[_]^n_->0,d[_]^n_->0,e[_]^n_->0,f[_]^n_->0,a[l1_]*a[l2_]->0,b[l1_]*b[l2_]->0,c[l1_]*c[l2_]->0,d[l1_]*d[l2_]->0,e[l1_]*e[l2_]->0,f[l1_]*f[l2_]->0};
        FullState=Expand[(DoublePaths/2)*v[0]];
      
        (*Now condition on four-folds, getting all possible terms*)
        TriggerableState=Expand[FullState*ZERO]/.{ZERO*a[dda_]*b[ddb_]*c[ddc_]*d[ddd_]*v[0]->a[dda]*b[ddb]*c[ddc]*d[ddd]*v[0]}/.{ZERO->0};
        AllABCCombinations=DeleteDuplicates[TriggerableState/.{Plus->List}/.{ww[_]->1}];
        NormCombTmp=AllABCCombinations/.{a[_]->1,b[_]->1,c[_]->1,d[_]->1,v[_]->1};
        NormComb=ConstantArray[1,Length[AllABCCombinations]]/NormCombTmp;
        AllABCCombinations=AllABCCombinations*NormComb;
      
        (*Creating all possible paths for this SRV*)
        AllEquations={};
        TargetEquations={};
        ExpandedStateZERO=Expand[TriggerableState*ZERO];
        For[ii=1,ii<=Length[AllABCCombinations],ii++,
          newEq=ExpandedStateZERO/.{ZERO*AllABCCombinations[[ii]]->1}/.{ZERO->0};
            For[jj=1,jj<=Length[AllStates],jj++,
              If[(AllABCCombinations[[ii]]==a[0]*b[AllStates[[jj,1]]]*c[AllStates[[jj,2]]]*d[AllStates[[jj,3]]]*v[0]),(*This term is part of the quantum state we want*)
                AppendTo[TargetEquations,newEq];
              ];
            ];
            AppendTo[AllEquations,newEq];
          ];
          AllEquations=AllEquations/.{a[_]->0,b[_]->0,c[_]->0,d[_]->0,e[_]->0,f[_]->0,g[_]->0,h[_]->0,i[_]->0,j[_]->0,a[_]^n_->0,b[_]^n_->0,c[_]^n_->0,d[_]^n_->0,e[_]^n_->0,f[_]^n_->0,g[_]^n_->0,h[_]^n_->0,i[_]^n_->0,j[_]^n_->0};
          TargetEquations=TargetEquations/.{a[_]->0,b[_]->0,c[_]->0,d[_]->0,e[_]->0,f[_]->0,g[_]->0,h[_]->0,i[_]->0,j[_]->0,a[_]^n_->0,b[_]^n_->0,c[_]^n_->0,d[_]^n_->0,e[_]^n_->0,f[_]^n_->0,g[_]^n_->0,h[_]^n_->0,i[_]^n_->0,j[_]^n_->0};

          (* Run the Optimization *)
          NormalisationConstant=Total[Abs[AllEquations]^2];
          Fidelity=Total[Abs[TargetEquations]]^2/(Length[TargetEquations]*NormalisationConstant);(*Definition of Fidelity*)
          vars=Sort[DeleteDuplicates[Cases[NormalisationConstant,_ww,{0,Infinity}]]];
          LenOfVariables=Length[vars];

          alpha=0.0; (*no L1 regularisation when we are just discovering states. Only important when using topological optimization*)
          Loss1=(1-Fidelity);
          Loss2=(1-Fidelity)+alpha*Total[Abs[vars]];

          fullvars=Transpose[{vars,RandomReal[{-1,1},Length[vars]]}];(*initial random values of weights*)
          sol=FindMinimum[Loss2,fullvars,AccuracyGoal->3,PrecisionGoal->3];(* Thank you Broyden, Fletcher, Goldfarb, Shanno*)

          L1val=Loss1/.sol[[2]];(*pure (1-fidelity) withough L1 term*)
          saveCount++;
        ];
        EndTime=AbsoluteTime[];

        If[saveCount<=MaxSaveCount,
          PrintLog[ToString[localDim[[2;;4]]]<>": "<>ToString[L1val,InputForm]<>" ("<>ToString[N[EndTime-StartTime],InputForm]<>" sec.), "<>ToString[saveCount]];
          NumOfGoodStates++;
          ,
          PrintLog["Final - "<>ToString[localDim[[2;;4]]]<>": No State ("<>ToString[N[EndTime-StartTime],InputForm]<>" sec.)"];
          NumOfBadStates++;
        ];
      ];
    ];
  ];
];

FinalTime=N[AbsoluteTime[]-TotalTimeStart];

PrintLog["\nNumOfGoodStates: "<>ToString[NumOfGoodStates]];
PrintLog["NumOfBadStates: "<>ToString[NumOfBadStates]];
PrintLog["FinalTime: "<>ToString[FinalTime]];
