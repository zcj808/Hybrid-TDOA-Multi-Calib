function [arrs_tdoas] = get_tdoa(file_id,thres,label)
arrs_tdoas=[];
arr_M=[6,6,6];
if label=="my"
    arrs_tdoas=[];
    for ii=1:length(arr_M)
        fs=16000;
        [b,a] = butter(6,[500,4000]/(fs/2));
        arr_tdoas=[];
        for i=1:arr_M(ii)
            id=i;
            filename = sprintf('audio/%d/arr%d-%d.wav',file_id,ii,id);
            [y,fs] = audioread(filename);
            f_data = filter(b,a,y');
            wid=100;
            step=50;
            scale_chirp_len=fix(0.6*fs);
            frame_num=fix((size(f_data,2)-wid)/step+1);
            rough_id=0;
            k=1;
            arrs_tdoa=[];
            energys=[];
            while k<frame_num-2
                energy=sum(f_data(step*(k-1)+1:step*(k-1)+wid).^2);
                energys=[energys,energy];
                if k>1 && energy-last_energy>thres && last_energy>0
                    rough_id=step*(k-1)+1-3*wid;
                    k=k+fix(scale_chirp_len/step);
                    if last_rough_id>0 && rough_id+scale_chirp_len<size(f_data,2)
                        ref_data=f_data(last_rough_id:last_rough_id+scale_chirp_len);
                        curr_data=f_data(rough_id:rough_id+scale_chirp_len);
                        tau = gccphat(curr_data',ref_data',fs);
                        X=fft(curr_data);
                        Y=fft(ref_data);
                        G=X.*conj(Y); % GCC-PHAT
                        gcc=G./abs(G);
                        gcc=ifftshift(ifft(gcc));
                        [max_v,max_id]=max(abs(gcc));
                        arrs_tdoa=[arrs_tdoa,(max_id+rough_id-last_rough_id-scale_chirp_len/2)/fs];
                    end
                else
                    k=k+1;
                end
                last_energy=energy;
                last_rough_id=rough_id;
            end
            sprintf("tdoa-s %d %d %d",file_id,ii,id)
            plot(energys)
            size(arrs_tdoa,2)
            arr_tdoas=[arr_tdoas;arrs_tdoa];
        end
        arrs_tdoas=[arrs_tdoas;median(arr_tdoas)];
    end
elseif label=="su" % select one mic. in each mic. array
    for ii=1:length(arr_M)
        fs=16000;
        [b,a] = butter(6,[500,4000]/(fs/2));
        arr_ids=[];
        f_datas=[];
        datas=[];
        for i=1:arr_M(ii)
            id=i;
            filename = sprintf('audio/%d/arr%d-%d.wav',file_id,ii,id);
            [y,fs] = audioread(filename);
            f_data = filter(b,a,y');
            f_datas = [f_datas;f_data];
            datas=[datas;y'];
            wid=100;
            step=50;
            scale_chirp_len=fix(0.6*fs);
            frame_num=fix((size(f_data,2)-wid)/step+1);
            last_energy=0;
            j=1;
            arr_id=[];
            energys=[];
            while j<frame_num-2
                energy=sum(f_data(step*(j-1)+1:step*(j-1)+wid).^2);
                if j>1 && energy-last_energy>thres && last_energy>0
                    rough_id=step*(j-1)+1-3*wid;
                    j=j+fix(scale_chirp_len/step);
                    arr_id=[arr_id,rough_id];
                else
                    j=j+1;
                end
                energys=[energys;energy];
                last_energy=energy;
            end
            sprintf("tdoa-m %d %d %d",file_id,ii,id)
            size(arr_id)
            arr_ids=[arr_ids;arr_id];
        end
        arr_data(ii).data=datas;
        arr_f_data(ii).data=f_datas;
        arr_arr_id(ii).id=arr_ids;
    end
    save(sprintf("data/%d/doa_audio.mat",file_id),"arr_data");
    save(sprintf("data/%d/doa_audio_id.mat",file_id),"arr_arr_id");
    su_tdoa_num=size(arr_arr_id(1).id,2);
    arrs_tdoas=[];
    for ii=2:length(arr_M)
        arr_tdoas=[];
        for i=1:arr_M(1)
            pat_datas=arr_f_data(1).data(i,:);
            pat_ids=arr_arr_id(1).id(i,:);
            for j=1:arr_M(ii)
                rev_datas=arr_f_data(ii).data(j,:);
                rev_ids=arr_arr_id(ii).id(j,:);
                tdoa=[];
                for kk=1:su_tdoa_num
                    pat_id=pat_ids(kk);
                    pat_data=pat_datas(pat_id:pat_id+scale_chirp_len);
                    rev_id=rev_ids(kk);
                    rev_data=rev_datas(rev_id:rev_id+scale_chirp_len);
                    tau = gccphat(rev_data',pat_data',fs);
                    X=fft(rev_data);
                    Y=fft(pat_data);
                    G=X.*conj(Y); % GCC-PHAT
                    gcc=G./abs(G);
                    gcc=ifftshift(ifft(gcc));
                    [max_v,max_id]=max(abs(gcc));
                    tdoa=[tdoa,tau+(rev_id-pat_id)/fs];
                end
                arr_tdoas=[arr_tdoas;tdoa];
            end
        end
        arr_tdoas=median(arr_tdoas);
        arrs_tdoas=[arrs_tdoas;arr_tdoas];
    end
end