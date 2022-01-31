a=`head -6 DOSCAR|tail -1|awk '{print $3}'` # 从DOSCAR中获取NEDOS的值
b=$((a + 6))    # 确定最后的行数
f=`awk '{if(NR==6)print $4}' DOSCAR`  # 从DOSCAR中获取费米能
sed -n '7,'$b' p' DOSCAR > DOS.dat  # 获取dos数据
awk '{print $1-'$f',$2}' DOS.dat > DOS-final.dat  # 减去费米能