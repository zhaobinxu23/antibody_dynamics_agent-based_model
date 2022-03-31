function MM = Trace_contact(M,Number_list)
for i = 1:length(M)
    if M(i,Number_list)> 0
        M(i,:) = zeros(1,length(M));
    end
end
 M(Number_list,:) = zeros(1,length(M));
MM = M;
end

    