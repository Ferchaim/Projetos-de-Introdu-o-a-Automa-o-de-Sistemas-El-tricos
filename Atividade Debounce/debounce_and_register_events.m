function [debounced_data, event_log] = debounce_and_register_events(data, duration, sampling_period)
    % Inicializacao
    debounced_data = zeros(size(data));
    buffer = zeros(size(data));
    counter = 0;
    event_log = [];

    % Loop pelos dados
    for i = 1:numel(data)
        if data(i) ~= buffer(i)
            counter = counter + 1;
        else
            counter = 0;
        end

        if counter >= duration
            debounced_data(i) = data(i);
        elseif counter == 0
            debounced_data(i) = data(i);
        else
            debounced_data(i) = debounced_data(i-1);
        end

        buffer(i) = data(i);

        % Registro de eventos
        if i > 1 && debounced_data(i) ~= debounced_data(i-1)
            event_type = debounced_data(i);
            event_time = (i - 1) * sampling_period;
            event_log = [event_log; event_time, event_type];
        end
    end
end
