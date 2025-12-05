function result = create_studies_table_burdenEM(...
    titles, prevalences, save_path, model_path, num_replicates, h2BurdenTrue)
%create_studies_table_burdenEM generates a 'studies' table for burdenEM
%   This table contains one row for each simulation that was performed

ii = 1;
for rep = 1:num_replicates
    for titles_idx = 1:length(titles)
        simulation = [titles{titles_idx},'.rep=',num2str(rep)];
        result.identifier{ii,1} = simulation;
        result.abbreviation{ii,1} = titles{titles_idx};
        if prevalences(titles_idx) == 0
            result.trait_type{ii,1} = 'continuous';
        else
            result.trait_type{ii,1} = 'binary';
        end
        result.description{ii,1} = '';
        result.sumstats_filename_pattern{ii,1} = ...
            [save_path, simulation, '.txt'];
        result.model_filename{ii,1} = ...
            [model_path, simulation, '.rds'];
        result.dataset{ii,1} = 'simulated';
        ii = ii+1;
    end
end


end