function Population = EnvironmentalSelection(Population,N,k)

    while length(Population) > N
        %% Diversity_calculation
        dist = pdist2(Population.decs,Population.decs);
        dist = sort(dist);
        for i=1:length(Population)
            d_var(i,:) = sum(dist(i,1:k));
        end

        dist = pdist2(Population.objs,Population.objs);
        dist = sort(dist);
        for i=1:length(Population)
            d_obj(i,:) = sum(dist(i,1:k));
        end

        m_var = sum(d_var)/length(Population);
        m_obj = sum(d_obj)/length(Population);

        for i=1:length(Population)
            f_var(i,:) = 1/(1+sum(d_var(i,:))/m_var);
            f_obj(i,:) = 1/(1+sum(d_obj(i,:))/m_obj);
            f_var_obj(i,:) = 1/(1+sum(d_var(i,:))/m_var+sum(d_obj(i,:))/m_obj);
        end

        z=min(f_var,[],1);
        Z=max(f_var,[],1);
        n_f_var=(f_var-z)./repmat(Z-z,length(f_var),1);

        z=min(f_obj,[],1);
        Z=max(f_obj,[],1);
        n_f_obj=(f_obj-z)./repmat(Z-z,length(f_obj),1);

        z=min(f_var_obj,[],1);
        Z=max(f_var_obj,[],1);
        n_f_var_obj=(f_var_obj-z)./repmat(Z-z,length(f_var_obj),1);
        
        Pd = zeros(1,length(Population));
        for i = 1:length(Population)
            if n_f_var > mean(n_f_var)
                Pd(i) = max(n_f_var(i), n_f_obj(i));
            elseif n_f_obj > mean(n_f_obj)
                Pd(i) = max(n_f_var(i), n_f_obj(i));
            elseif n_f_var < mean(n_f_var)
                Pd(i) = min(n_f_var(i), n_f_obj(i));
            elseif n_f_obj < mean(n_f_obj)
                Pd(i) = min(n_f_var(i), n_f_obj(i));
            else
                Pd(i) = n_f_var_obj(i);
            end
        end

        [~,x] = min(Pd);
        Population(x) = [];
    end