function [MFR_data, Isp_data, PM_data, DM_data, H_data, MV_data, TB_data, F_data, IT_data] = createData(A)
    index = 1;
    for i = 1:length(A(:,1,1,1,1))
        for j = 1:length(A(1,:,1,1,1))
            for k = 1:length(A(1,1,:,1,1))
                for l = 1:length(A(1,1,1,:,1))
                    MFR_data(index) = A(i,j,k,l,1);
                    Isp_data(index) = A(i,j,k,l,2);
                    PM_data(index) = A(i,j,k,l,3);
                    DM_data(index) = A(i,j,k,l,4);
                    H_data(index) = A(i,j,k,l,5);
                    MV_data(index) = A(i,j,k,l,6);
                    TB_data(index) = A(i,j,k,l,7);
                    F_data(index) = A(i,j,k,l,8);
                    IT_data(index) = A(i,j,k,l,9);
                    index = index + 1;
                end
            end
        end
    end
end