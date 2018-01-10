using MendelLocationScores
using MendelBase
using Search
using SearchSetup
using DataFrames

# @testset "initialize_optimization_two_point_linkage" begin

# end

@testset "prior function" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["flanking_distance"] = [0.5, 0.5]
  	keyword["flanking_markers"] = 1
    keyword["gender_neutral"] = true
    keyword["lod_score_table"] = "Lod_Score_Frame.txt"

    keyword["parameters"] = 1
    keyword["points"] = 10
    keyword["eliminate_genotypes"] = true
  	keyword["lump_alleles"] = true
    parameter = set_parameter_defaults(keyword)
    process_keywords!(keyword, "location scores Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)
    keyword["eliminate_genotypes"] = true
    keyword["lump_alleles"] = true

    par = parameter.par #this is [0.0]
    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree,
        person, nuclear_family, locus, keyword)

    # use dictionaries to assign prior probabilities
    # the numbers are provided in the locus frame
    dic_locus1 = Dict("1" => 0.62, "2" => 0.17, "3" => 0.14, "4" => 0.07)
    dic_locus2 = Dict("1" => 0.003, "2" => 0.997)
    locus_3_sum = 0.41 + 0.14 + 0.03 + 0.39 #need to normalize allele freq to 1 at locus 3
    dic_locus3 = Dict("1" => 0.41/locus_3_sum, "2" => 0.14/locus_3_sum,
        "3" => 0.03/locus_3_sum, "4" => 0.39/locus_3_sum)

    # loop through to compute probabilities
    # n is the variable iterating from instruction.start[ped] to instruction.finish[ped].
    for ped in 1:pedigree.pedigrees
        for n = instruction.start[ped]:instruction.finish[ped]-1
            operation = instruction.operation[n]
            start = instruction.extra[n][1]
            finish = instruction.extra[n][2]
            i = instruction.extra[n][3]
            if operation != penetrance_and_prior_array continue end #avoids some array access errors
            if person.mother[i] != 0 continue end #prior prob doesnt exist for non founder

            #
            # Construct the parent's multiple locus genotypes.
            #
            genotypes = MendelBase.genotype_count(person, locus, i, start, finish)
            multi_genotype = MendelBase.construct_multigenotypes(person, locus, start, finish,
                                                genotypes, i)

            for j = 1:genotypes
                prob = MendelLocationScores.prior_location_score(person, locus, 
                    multi_genotype[:, :, j], par, keyword, start, finish, i)
                answer = 1.0

                # tally probabilities contributed by the 3 locus 
                answer *= dic_locus1[string(multi_genotype[1, 1, j])]
                answer *= dic_locus1[string(multi_genotype[2, 1, j])]

                answer *= dic_locus2[string(multi_genotype[1, 2, j])]
                answer *= dic_locus2[string(multi_genotype[2, 2, j])]

                answer *= dic_locus3[string(multi_genotype[1, 3, j])]
                answer *= dic_locus3[string(multi_genotype[2, 3, j])]

                @test answer == prob
            end
        end
    end
end

@testset "penetrance function" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["flanking_distance"] = [0.5, 0.5]
  	keyword["flanking_markers"] = 1
    keyword["gender_neutral"] = true
    keyword["lod_score_table"] = "Lod_Score_Frame.txt"

    keyword["parameters"] = 1
    keyword["points"] = 10
    keyword["eliminate_genotypes"] = true
  	keyword["lump_alleles"] = true
    parameter = set_parameter_defaults(keyword)
    process_keywords!(keyword, "location scores Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)
    keyword["eliminate_genotypes"] = true
    keyword["lump_alleles"] = true

    par = parameter.par #this is [0.0]
    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree,
        person, nuclear_family, locus, keyword)

    # n is the variable iterating from instruction.start[ped] to instruction.finish[ped].
    for ped in 1:pedigree.pedigrees
        for n = instruction.start[ped]:instruction.finish[ped]-1
            operation = instruction.operation[n]
            start = instruction.extra[n][1]
            finish = instruction.extra[n][2]
            i = instruction.extra[n][3]
            if operation != penetrance_and_prior_array continue end #this avoids some array access errors

            #
            # Construct the parent's multiple locus genotypes.
            #
            genotypes = MendelBase.genotype_count(person, locus, i, start, finish)
            multi_genotype = MendelBase.construct_multigenotypes(person, locus, start, finish,
                                                genotypes, i)

            for j = 1:genotypes
                number = MendelLocationScores.penetrance_location_score(person, locus, 
                    multi_genotype[:, :, j], par, keyword, start, finish, i)
                @test number == 1.0 #this is always 1 
            end
        end
    end
end

# @testset "transmission function" begin

# end

@testset "wrapper and basics" begin

end


