# 平面三角形单元
# 单位：m,N
wipe
puts "System"
# 二维，每个节点2个自由度
model basic -ndm 2 -ndf 2
puts "node"

node 1 2.0 1.0
node 2 2.0 0.0
node 3 0.0 1.0
node 4 0.0 0.0

puts "restraint"
fix 3 1 1;
fix 4 1 1;

puts "material"
# 多轴材料，弹性模量，泊松比。
set E 1e7 
set NU 0.33333
nDMaterial ElasticIsotropic 1 $E $NU

puts "element"
# 桁架单元：单元编号，3个节点（逆时针），type( "PlaneStrain" or "PlaneStress")，材料编号
set thickness 0.1
element tri31 1 2 3 4 $thickness PlaneStress 1
element tri31 2 2 1 3 $thickness PlaneStress 1

puts "recorder"
# 输出1-4号节点的位移
recorder Node -file node0.out -time -nodeRange 1 4 -dof 1 2 disp


puts "loading"

pattern Plain 1 Linear {
	# 点荷载， FX，FY
load 1 0.0 -50000.0 
load 2 0.0 -50000.0 
}

puts "analysis"
constraints Plain
numberer Plain
system BandGeneral
test EnergyIncr 1.0e-6 200
algorithm Linear
integrator LoadControl 1
analysis Static
analyze 1
