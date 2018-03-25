#writedlm("../Single_cell_data/sample_data_summary.txt", df_stats) 

#=
open("../Single_cell_data/sample_data_summary.txt") do df_stats
    for i in enumerate(eachline(df_stats))
      println(i)
    end
end
=#

f = open("../Single_cell_data/sample_data_summary.txt");
for ln in eachline(f)
       print("$(length(df_stats)), $ln")
end
close(f)
