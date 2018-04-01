function Triangle2D3Node_Display
global gNode gElement gMaterial gBC1 gNF gDF gK gDelta
TR=triangulation(gElement(:,1:3),gNode);
triplot(TR);
return