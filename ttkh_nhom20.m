clear
clc

thoat = 0;
while (thoat == 0)
% tao menu lua chon cac kieu do thi.
s = menu('MO PHONG SONG NUOC','Mo hinh theo ham sin','Giao thoa','Co bien','Co vach ngan','Mua roi','thoat'); 
switch s
    case 1
        x = 0:.5:100;   
        y = 0:.5:100;
        [X,Y] = meshgrid(x,y); % tao mot luoi diem, xem moi phan tu trong luoi la ot phan tu nuoc
        z = [];
        z(1:length(x),1:length(y)) = 0; %tao ma tran z chua toa do cua cac phan tu song
        n = 0; 
        t = 0; 
        while n < 300  
            for i = 1:length(x)
                for j = 1:length(y)
                    d = sqrt((x(i)-50)^2+(y(j)-50)^2);%d la khoang cach tu nguon den phan tu dang xet
                    if (d>t) % d > v*t voi v = 1(m/s), v*t la ban kinh song truyen toi trong thoi gian t
                        z(i,j) = 0; % neu song chua truyen toi thi phan tu do chua dao dong
                    else
                        z(i,j) = 0.5*cos(2*pi*t/3+2*pi*d/3); 
                        % khi song da truyen toi thi phan tu nuoc do dao
                        % dong voi pt: x = 0.5*cos(2*pi*t/3+2*pid/3).
                        % lay T = 1/3(s), lamda = 1/3(m).
                    end
                end
            end
        pause(0.000000001);% thoi gian dung man hinh 
        t = t + 0.1; % thoi gia tang 0.1s sau moi vong lap
        surfl(z); % ve mat cua ma tran z
        zlim([-5,5]); % truc z trong do thi
        view(30,80); %goc nhin 
        colormap winter; % tao mau winter cho do thi
        shading interp; % kieu do thi
        n = n+1; 
        end
    case {2,3,4,5}
        nx = 251;
        ny = 251; 
        k = 2;
        [x y] = meshgrid(linspace(0,250,nx),linspace(0,250,ny));
        % tao mot luoi diem, xem moi phan tu trong luoi la ot phan tu nuoc
        % ham linspace tao mot vecto co 250 phan tu trong khoang 0->nx
        wave = zeros(nx,ny); % tao ma tran 0, bieu dien cac phan tu song
        p = zeros(nx,ny); % tao ma tran 0
        q = zeros(nx,ny); % tao ma tran 0
        wall = zeros(nx,ny); % tao ma tran 0, bieu dien cac vach ngan
        switch s
            case 2
                source_number = 2; % so nguon
                source = [100 70; 800 100];% toa do cua nguon
                border_number = 4;% so vach ngan
                border = [0 1 1 251; 1 1 1 251; 0 251 1 251; 1 251 1 251];
            case 3
                source_number = 1;%so nguon
                source = [125 125]; % toa do cua nguon
                border_number = 4;% so vach ngan
                border = [0 1 1 251; 1 1 1 251; 0 251 1 251; 1 251 1 251];
            case 4
                source_number = 2;% soo nguon
                source = [150 50; 100 50];% toa do cua nguon
                border_number = 6;% so vach ngan
                border = [0 1 1 251;1 1 1 251;0 251 1 251;1 251 1 251;0 100 1 110;0 100 140 250];
            case 5
                source_number = 0; % so nguon
                border_number = 4; % so vach ngan
                border = [0 1 1 251; 1 1 1 251; 0 251 1 251; 1 251 1 251];
        end
        
        for h = 1:source_number % lap source_number lan voi source_number la so nguon
            for i=1:nx 
                for j=1:ny
                    wave(i,j) = wave(i,j) + 10*exp((-((i-source(h,1))^2 + (j-source(h,2))^2))/(k^2));
                    % tao ma tran song voi cac nguon 
                end
            end
        end
        for h = 1:border_number % lap border_number lan voi border_number la so vach ngan
            % tao cac vach ngan co do cao bang 2
            if (border(h,1) == 0) 
                wall(border(h,3):border(h,4),border(h,2)) = 2;
            else
                wall(border(h,2),border(h,3):border(h,4)) = 2;
            end
        end
        surf(x,y,wave);% ve mat 
        view(30,60); % goc nhin
        
        p = wave;
        n = 0;
        dt=0.3;
        rt = 0;
        v = 1;
        while n < 300
            k = 2;
            if (rt<2e-2)
               rt = rt + dt*(3e-5);
               v = exp(-rt); % v = e^(-rt)
            end

            for i=2:nx-1
                for j=2:ny-1
                    if (wall(i,j)) 
                        continue;
                    end
                    q(i,j)=(2*wave(i,j)-p(i,j)+0.5*(wave(i+1,j)+wave(i-1,j)+wave(i,j+1)+wave(i,j-1)-4*wave(i,j)))*v;
                    % tinh toa do cua cac phan tu song
                end
            end
            p = wave;
            wave = q;
    
            if (s == 5 && mod(n,30) == 0) % mua roi, sau 30 vong lap se co 1 giot mua
                k = 4;
                while 1
                    x0 = randi([1 nx],1,1);
                    y0 = randi([1 ny],1,1);
                    %tao random toa do (x0,y0) la vi tri mua roi
                    if (not(wall(x0,y0))) 
                        break;
                    end
                end
                for i=1:nx 
                    for j=1:ny
                        if (wall(i,j))
                        eta(i,j) = 0;
                        else
                            p(i,j) = p(i,j) - 5*exp((-((i-x0)^2 + (j-y0)^2))/(k^2));
                            wave(i,j) = wave(i,j) - 5*exp((-((i-x0)^2 + (j-y0)^2))/(k^2));
                            % tinh toa do cua cac phan tu song
                        end
                    end
                end
            end
            figure(1); %T?o m?i hình ?nh (?? th?).
            F = wave + wall;
            surf(x,y,F) % ve mat
            zlim([-5 5]); % truc z trong do thi
            view(30,60); % goc nhin
            if s == 5 
                colormap winter; % tao mau winter cho do thi mua roi
            else
                colormap Gray;  % tao mau Gray cho do thi con lai
            end
            shading interp % kieu do thi
            n = n + 1;
        end

    case 6
        thoat = 1;
end
end
    