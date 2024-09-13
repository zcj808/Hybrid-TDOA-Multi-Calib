%%
for i=1:15
    for j=1:3
        for k=1:6
            filename=sprintf('audio/%d/arr%d-%d.wav',i,j,k);
            filename
            [y,fs] = audioread(filename);
%             if k==1
%                 filename = sprintf('audio/%d/arr%d-6.wav',i,j);
%             else
%                 filename = sprintf('audio/%d/arr%d-%d.wav',i,j,k-1);
%             end
%             filename
            audiowrite(filename,y,fs,'BitsPerSample',16);
        end
    end
end
%%
clear
clc
thres=1;
for i=1:15
    i
    for j=1:3
        pcm_file = sprintf('audio3/%d/mic_demo_vvui_ori%d.pcm',i,j);
        % read .pcm file
        M=6; % mic. number
        fid = fopen(pcm_file, 'rb');
        data = fread(fid, [M+2,Inf], 'int32'); % determine whether decoding is correct
        fclose(fid);
        gain=2^28;
        data=double(single(data)./gain); % data(7,:) and data(8,:) are reference channels
        fs=16000;
        [b,a] = butter(6,[500,4000]/(fs/2));
        rough_id=0;
        for kk=1:6
            f_data = filter(b,a,data(kk,:));
            plot(f_data);
            wid=100;
            step=50; % windows length equals to step length
            scale_chirp_len=fix(0.6*fs); % unit of 0.5 is second
            chirp_len=fix(0.3*fs);
            frame_num=fix((size(f_data,2)-wid)/step+1);
            if kk==1
                k=1;
                arrs_tdoa=[];
                energys=[];
                while k<frame_num-2
                    energy=sum(f_data(step*(k-1)+1:step*(k-1)+wid).^2);
                    energys=[energys,energy];
                    if k>1 && energy-last_energy>thres && last_energy>0
                        rough_id=step*(k-1)+1-3*wid-floor(1000*rand)-5*fs;
                        break
                    else
                        k=k+1;
                    end
                    last_energy=energy;
                end
            end
            y=data(kk,rough_id:end);
            plot(y);
            if kk==1
                filename = sprintf('audio/%d/arr%d-6.wav',i,j);
            else
                filename = sprintf('audio/%d/arr%d-%d.wav',i,j,kk-1);
            end
%             filename
%             filename = sprintf('audio/%d/arr%d-%d.wav',i,j,kk);
            audiowrite(filename,y,16000,'BitsPerSample',16);
        end
    end
end