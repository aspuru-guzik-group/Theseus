(*
    Conceptual understanding through efficient inverse-design of quantum optical experiments
    
    Example: two-qubit CNOT
    Mario Krenn (mario.krenn@utoronto.ca), Jakob Kottmann, Nora Tischler, Alán Aspuru-Guzik (alan@aspuru.com)
    https://arxiv.org/abs/2005.06443
    mariokrenn.wordpress.com/ https://www.matter.toronto.edu/
    
    We show how an experiment for a heralded CNOT can be translated to graphs.
    The experiment has been performed in 2004 (Phys. Rev. Lett. 93, 020504, 2004)
    by Gasparoni et.al. After translating the experiment to the graph, we plot the
    corresponding graph (requires Mathematica 12.1).
*)


(*
    Define optical elements that are used,
    Polarising Beam Splitter, Half-wave plate,
    polarising BS in the D/A basis,
    transformation between computation and D/A basis
*)

PBS[expr_,a_,b_]:=expr/.{a[0]->b[0],a[1]->a[1],b[0]->a[0],b[1]->b[1]}
HWP[expr_,a_]:=expr/.{a[0]->a[0]+a[1],a[1]->a[0]-a[1]}

PBSDA[expr_,a_,b_]:=expr/.{a[AA]->b[AA],a[DD]->a[DD],b[AA]->a[AA],b[DD]->b[DD]}

ToAD[expr_,a_]:=expr/.{a[0]->a[DD]+a[AA],a[1]->a[DD]-a[AA]}
ToHV[expr_,a_]:=expr/.{a[DD]->a[0]+a[1],a[AA]->a[0]-a[1]}


(*These three lists will contain the edge,color,weight information of the four generated Graphs*)
GraphEdges={};
GraphColors={};
GraphWeights={};

(*For each of the logical inputs 00,01,10,11, we calculate the setup and create a Graph*)
For[CC=0,CC<=1,CC++,
  For[TT=0,TT<=1,TT++,
    psi=Va[CC]*a[CC]+(b[0]*c[1]+b[1]*c[0])+Vd[TT]*d[TT];
    
    (*Setup from Gasparoni et.al. *)
    psi=ToAD[psi,c];
    psi=ToAD[psi,d];
    psi=PBSDA[psi,c,d];
    psi=ToHV[psi,c];
    psi=ToHV[psi,d];
    psi=PBS[psi,a,b];
    CurrentState=Expand[psi]/.{Plus->List};
    CurrentGraphNorm=CurrentState/.{x_[_]->1};
    AppendTo[GraphWeights,CurrentGraphNorm];(*edge weights*)
    
    CurrentStateNorm=CurrentState/CurrentGraphNorm;
    AppendTo[GraphEdges,CurrentStateNorm/.{x_[l1_]*y_[l2_]->{x,y}}];(*Extract the edges from states*)
    AppendTo[GraphColors,CurrentStateNorm/.{x_[l1_]*y_[l2_]->{l1,l2}}/.{0->Blue,1->Red,2->Green}];(*Get the colors from state*)
    
    psi=psi/.{a[0]->a[DD]+a[AA],a[1]->a[DD]-a[AA]};
    psi2=Expand[psi^3];
    
    psiHeralded=Expand[ZERO*psi2]/.{ZERO*Va[CC]*Vd[TT]*a[DD]*d[1]->1}/.{ZERO->0};
    Expand[psiHeralded];
    Print["|"<>ToString[CC]<>","<>ToString[TT]<>">   ->    "<>ToString[Expand[psiHeralded/12]]];(*we certify that CNOT works correctly*)
  ];
];

(*
    The following code plots the Graph.
    It requires Mathematica 12.1, because it uses EdgeTaggedGraph
    Code adapted from kglr at https://mathematica.stackexchange.com/a/222043/12750
*)

StoreGraph=Array["object",{Length[GraphEdges]}];

For[ii=1,ii<=Length[GraphEdges],ii++,
  edges=DirectedEdge@@@GraphEdges[[ii]];
  edgecolors=GraphColors[[ii]];
  
  taggededges=EdgeList@EdgeTaggedGraph@edges;
  coloring=AssociationThread[taggededges,edgecolors];
  
  eShapeFunction2=Module[{c=coloring@#2,bsf=BSplineFunction@#,s=Partition[Subdivide[Length@coloring@#2],2,1]},{CapForm["Butt"],Thread[{c,Line/@(bsf/@Subdivide[##,100]&@@@s)}]}]&;
  
  graphs=Graph[Reverse@{Va, Vd,a,b,c,d},taggededges,GraphLayout->{"CircularEmbedding","Offset"->2.5*(2*Pi/6)},VertexLabels->Placed["Name",Center],VertexSize->.4,VertexStyle->White,VertexLabelStyle->Directive[FontFamily->"Times",Large],EdgeStyle->Directive[CapForm["Round"],Opacity[.7],AbsoluteThickness[5]],PerformanceGoal->"Quality",EdgeShapeFunction->eShapeFunction2];
  StoreGraph[[ii]]=Show[graphs];
];
StoreGraph

