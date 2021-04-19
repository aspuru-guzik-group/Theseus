(*
    Conceptual understanding through efficient inverse-design of quantum optical experiments
    Example: Translating a Graph to Path encoding and  Bulk Optics
    Mario Krenn (mario.krenn@univie.ac.at), Jakob Kottmann, Nora Tischler, Alán Aspuru-Guzik (alan@aspuru.com)
    https://arxiv.org/abs/2005.06443
    mariokrenn.wordpress.com/ https://www.matter.toronto.edu/
    Toronto, 19.04.2021
*)

(*
    Define optical elements that are used, Polarising Beam Splitter, 
    Half-wave plate, polarising BS in the D/A basis, transformation 
    between computation and D/A basis
*)

PBS[expr_, a_, b_] := expr /. {a[0] -> b[0], a[1] -> a[1], b[0] -> a[0], b[1] -> b[1]}
BS[expr_, a_, b_] := expr /. {a[nn_] -> (I*a[nn] + b[nn])/Sqrt[2], b[nn_] -> (a[nn] + I*b[nn])/Sqrt[2]}
HWP1[expr_, a_] := expr /. {a[0] -> a[0] + a[1], a[1] -> a[0] - a[1]}
HWP2[expr_, a_] := expr /. {a[0] -> a[1], a[1] -> a[0]}
SPDC[expr_, a_, na_, b_, nb_] := expr + a[na]*b[nb]
Absorb[expr_, a_] := expr /. {a[n_] -> 0}
RenamePath[expr_, a_, a1_] := expr /. {a1[nn_] -> a[nn]}

AlphabetToNum := {a -> 1, b -> 2, c -> 3, d -> 4, e -> 5, f -> 6};
NumToAlphabet := {1 -> a, 2 -> b, 3 -> c, 4 -> d, 5 -> e, 6 -> f};



(*
    This Graph will be translated to a setup
*)
graphEdges={a[1]*b[1],a[1]*f[1],a[1]*f[0],c[1]*f[1],c[0]*d[0],c[1]*d[1]};

NumOfVertices=LetterNumber[ToString[Last[Sort[Flatten[graphEdges/.{x_[l1_]*y_[l2_]->{x,y}}]]]]];
NumOfDimensions=Length[DeleteDuplicates[Flatten[graphEdges/.{x_[l1_]*y_[l2_]->{l1,l2}}]]];

setupList={};

(*
    First, we create a number of photon pair sources
*)
For[ii=1,ii<=Length[graphEdges],ii++,
    edgeinfo=graphEdges[[ii]]/.{a_[na_]*b_[nb_]->{a,na,b,nb}};
    Num1=edgeinfo[[1]]/.AlphabetToNum;
    Num2=edgeinfo[[3]]/.AlphabetToNum;
    AppendTo[setupList,"SPDC[XXX,"<>ToString[edgeinfo[[1]]]<>ToString[ii]<>","<>ToString[edgeinfo[[2]]]<>","<>ToString[edgeinfo[[3]]]<>ToString[ii]<>","<>ToString[edgeinfo[[4]]]<>"]"];
];

LossList={};


(*
    Now we combine different paths with beam splitters
    (and add an absorber in the unused output)
*)
graphEdgesTmp=graphEdges;
For[vv=1,vv<=NumOfVertices,vv++,
    For[dd=0,dd<NumOfDimensions,dd++,
        AlreadyUsedName="";
        For[ii=1,ii<=Length[graphEdges],ii++,
            edgeinfo=graphEdges[[ii]]/.{a_[na_]*b_[nb_]->{a,na,b,nb}};
            If[((edgeinfo[[1]]==vv/.NumToAlphabet)&&(edgeinfo[[2]]==dd)),
                If[StringLength[AlreadyUsedName]>0,
                    AppendTo[setupList,"BS[XXX,"<>AlreadyUsedName<>","<>ToString[edgeinfo[[1]]]<>ToString[ii]<>"]"];
                    AppendTo[setupList,"Absorb[XXX,"<>ToString[edgeinfo[[1]]]<>ToString[ii]<>"]"];
                    lpath=ToExpression[ToString[edgeinfo[[1]]]<>ToString[ii]];
                    AppendTo[LossList,lpath];
                    ,
                    AlreadyUsedName=ToString[edgeinfo[[1]]]<>ToString[ii];
                ];
            ];
            
            If[((edgeinfo[[3]]==vv/.NumToAlphabet)&&(edgeinfo[[4]]==dd)),
                If[StringLength[AlreadyUsedName]>0,
                    AppendTo[setupList,"BS[XXX,"<>AlreadyUsedName<>","<>ToString[edgeinfo[[3]]]<>ToString[ii]<>"]"];
                    AppendTo[setupList,"Absorb[XXX,"<>ToString[edgeinfo[[3]]]<>ToString[ii]<>"]"];
                    lpath=ToExpression[ToString[edgeinfo[[3]]]<>ToString[ii]];
                    AppendTo[LossList,lpath];
                    ,
                    AlreadyUsedName=ToString[edgeinfo[[3]]]<>ToString[ii];
                ];
            ];
        ];
    ];
];
Print["setupList for path encoding (as it is usually done in on-chip devices): ",setupList];

(*
    Finally, for bulk optics, we combine the remaining paths with polarizing beam splitters
*)
For[vv=1,vv<=NumOfVertices,vv++,
    Pos0=0;
    Pos1=0;
    For[ii=1,ii<=Length[graphEdges],ii++,
    
        edgeinfo=graphEdges[[ii]]/.{a_[na_]*b_[nb_]->{a,na,b,nb}};
        
        If[(((edgeinfo[[1]]==vv/.NumToAlphabet)&&(edgeinfo[[2]]==0))||((edgeinfo[[3]]==vv/.NumToAlphabet)&&(edgeinfo[[4]]==0)))&&Pos0==0,
            Pos0=ii;
        ];
        
        If[(((edgeinfo[[1]]==vv/.NumToAlphabet)&&(edgeinfo[[2]]==1))||((edgeinfo[[3]]==vv/.NumToAlphabet)&&(edgeinfo[[4]]==1)))&&Pos1==0,
            Pos1=ii;
        ];
    ];
    
    path0=ToExpression[ToString[vv/.NumToAlphabet]<>ToString[Pos0]];
    path1=ToExpression[ToString[vv/.NumToAlphabet]<>ToString[Pos1]];
    
    If[Pos0>0&&Pos1>0&&!MemberQ[LossList,path0]&&!MemberQ[LossList,path1],
        newComb="PBS[XXX,"<>ToString[path0]<>","<>ToString[path1]<>"]";
        AppendTo[setupList,newComb];
    ];
];

CurrLength=Length[setupList];
For[ii=1,ii<=Length[graphEdges],ii++,
    For[jj=1,jj<=NumOfVertices,jj++,
        currLett=ToString[jj/.NumToAlphabet];
        AppendTo[setupList,"RenamePath[XXX,"<>currLett<>","<>currLett<>ToString[ii]<>"]"];
    ];
];
Print["\n\nsetupList for polarisation encoding (as it is usually done in bulk optics): ",setupList[[1;;CurrLength]]];


XXX=0;
For[jj=1,jj<=Length[setupList],jj++,
    XXX=ToExpression[setupList[[jj]]];
];
Print["\n\nResulting state (Phases and Amplitude can directly by adjusted by the laser which pumps the initial SPDC crystals): ",XXX]

