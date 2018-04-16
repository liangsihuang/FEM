# 平面三角形单元
# 单位：mm,N
wipe
puts "System"
# 二维，每个节点2个自由度
model basic -ndm 2 -ndf 2
puts "node"

node 1 0.0 0.0
node 2 500.0 0.0
node 3 1000.0 0.0
node 4 1500.0 0.0
node 5 2000.0 0.0
node 6 0.0 500.0
node 7 500.0 500.0
node 8 1000.0 500.0
node 9 1500.0 500.0
node 10 2000.0 500.0
node 11 0.0 1000.0
node 12 500.0 1000.0
node 13 1000.0 1000.0
node 14 1500.0 1000.0
node 15 2000.0 1000.0

puts "restraint"
fix 1 1 1;
fix 6 1 1;
fix 11 1 1;

puts "material"
# 多轴材料，弹性模量，泊松比。
set E 2.1e5 
set NU 0.3
nDMaterial ElasticIsotropic 1 $E $NU

puts "element"
# 桁架单元：单元编号，3个节点（逆时针），type( "PlaneStrain" or "PlaneStress")，材料编号
set thickness 3
element tri31 1 1 2 6 $thickness PlaneStress 1
element tri31 2 2 7 6 $thickness PlaneStress 1
element tri31 3 2 3 7 $thickness PlaneStress 1
element tri31 4 3 8 7 $thickness PlaneStress 1
element tri31 5 3 4 8 $thickness PlaneStress 1
element tri31 6 4 9 8 $thickness PlaneStress 1
element tri31 7 4 5 9 $thickness PlaneStress 1
element tri31 8 5 10 9 $thickness PlaneStress 1
element tri31 9 6 7 11 $thickness PlaneStress 1
element tri31 10 7 12 11 $thickness PlaneStress 1
element tri31 11 7 8 12 $thickness PlaneStress 1
element tri31 12 8 13 12 $thickness PlaneStress 1
element tri31 13 8 9 13 $thickness PlaneStress 1
element tri31 14 9 14 13 $thickness PlaneStress 1
element tri31 15 9 10 14 $thickness PlaneStress 1
element tri31 16 10 15 14 $thickness PlaneStress 1

puts "recorder"
# 输出1-4号节点的位移
recorder Node -file node0.out -time -nodeRange 1 15 -dof 1 2 disp

recorder Element -file ele0.out -time -eleRange 1 16 stresses
puts "loading"

pattern Plain 1 Linear {
	# 点荷载， FX，FY
load 15 0.0 -1000.0 
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
