classdef rgasp_model
    %RGASP_MODEL RGASP model definition

    properties
        model
        samples
    end

    methods
        function obj = rgasp_model(model)
            obj.model = model;
        end

        function pred = predict(obj, x_new, options)
            arguments
                obj
                x_new
                options.idxSamples = nan;
                options.B = 500;
            end
            idxSamples = options.idxSamples;
            if isnan(idxSamples) 
                idxSamples = 1:options.B;
            end

            pred = rgasp_simulate(obj.model, x_new, options.B);
            pred = pred(:, idxSamples)';

        end
    end
end