
function out=groupingPreferenceMaker(N, K, ZipfParameter, group_count, group_interest)
    popularityArray = zeros(K, N);
    group_interest(1,3) = group_interest(3) / (group_count - 2);
    %group_interest
    group_size = (N / group_count);
    popularity = zeros(group_count, group_size);
    sum_popularity = zeros(group_count, group_size);
    allziph=zeros(1,group_count);
    for group_idx=1 : group_count
        for n=1 : group_size
            allziph(1,group_idx)=allziph(1,group_idx)+ ((1/n)^ZipfParameter(group_idx));
        end
    end
    for group_idx=1 : group_count
        %sumAllProbe = 0;
        for idx=1 : group_size
            popularity(group_idx,idx) = (1 / idx)^ZipfParameter(group_idx)/allziph(1,group_idx);
            %sumAllProbe = sumAllProbe + popularity(group_idx,idx);
            if idx > 1
                sum_popularity(group_idx,idx) = sum_popularity(group_idx,idx - 1) + popularity(group_idx,idx);
            else
                sum_popularity(group_idx,1) = popularity(group_idx,1);
            end
        end
        %sum_popularity(group_idx) = sum_popularity(group_idx) / sumAllProbe;
        %popularity(group_idx) = popularity(group_idx) / sumAllProbe;
    end
    %popularity
    for k=1 : K
%         Prob = zeros(1,N);
%         main_group = int32((group_count-1) * rand)+1;
%         main_group
%         slave_group = int32((group_count-1) * rand)+1;
%         if main_group == slave_group
%             slave_group =slave_group+1;
%             if slave_group > group_count
%                 slave_group=1;
%             end
%         end
%         slave_group
        Prob = zeros(1,N);
        FirstPopular_group = mod(k,group_count-1)+1;
        SecondPopular_group = int32((group_count-1) * rand)+1;
        if SecondPopular_group == FirstPopular_group
             SecondPopular_group =SecondPopular_group+1;
             if SecondPopular_group > group_count
                 SecondPopular_group=1;
             end
        end
        if SecondPopular_group > group_count
              SecondPopular_group=1;
        end
        idx = 1;
        for group_index=1 : group_count
            if group_index == FirstPopular_group
                for index=1 : group_size
                    Prob(1,idx) = popularity(group_index,index) * group_interest(1,1);
                    idx = idx+1;
                end
            elseif group_index == SecondPopular_group
                for index=1 : group_size
                    Prob(1,idx) = popularity(group_index,index) * group_interest(1,2);
                    idx = idx+1;
                end
            else
                for index=1 : group_size
                    Prob(1,idx) = popularity(group_index,index) * group_interest(1,3);
                    idx = idx+1;
                end
            end
        end
        popularityArray(k,:) = Prob(1,:);
    end
    %popularityArray
    out= popularityArray;
end
