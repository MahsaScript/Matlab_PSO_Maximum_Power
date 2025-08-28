% تنظیمات اولیه IDPSO
num_particles = 50;  % تعداد ذرات
max_iter = 4000;     % تعداد تکرارها
w_max = 0.9;         % حداکثر وزن اینرسی
w_min = 0.4;         % حداقل وزن اینرسی
c1 = 2.0;            % ضریب شناختی
c2 = 2.0;            % ضریب اجتماعی

% اطلاعات ورودی از برگه داده
Voc = 32.2;          % ولتاژ مدار باز (ولت)
Isc = 8.21;          % جریان اتصال کوتاه (آمپر)
VMPP = 26.3;         % ولتاژ در نقطه حداکثر توان (ولت)
IMPP = 7.61;         % جریان در نقطه حداکثر توان (آمپر)
Nc = 36;             % تعداد سلول‌های متصل سری

% محدوده جستجوی پارامترها
bounds = [0, 0.5;   % مقاومت سری (Rs)
          0, 1000;  % مقاومت موازی (Rp)
          1e-10, 1e-6]; % جریان اشباع معکوس دیود (Io)

% تابع هدف برای کاهش خطا
fitness_function = @(params) calculate_error(params, Voc, Isc, VMPP, IMPP, Nc);

% مقداردهی اولیه موقعیت‌ها و سرعت‌ها
% positions = bounds(:,1) + (bounds(:,2) - bounds(:,1)) .* rand(num_particles, 3);

positions = bounds(:,1) + (bounds(:,2) - bounds(:,1)).* rand() ;
% a = (bounds(2,:) - bounds(1,:));
% b = bounds(1,:) +a;
% c=rand(num_particles, 2)';
% positions =b * c;
% positions=positions';
% positions = bounds(1,:) + (bounds(2,:) - bounds(1,:)) .* rand(num_particles, 2)';

velocities = rand(num_particles, 3);

% مقداردهی اولیه بهترین موقعیت‌ها
p_best_positions = positions;
p_best_values = inf(num_particles, 1);
g_best_position = [];
g_best_value = inf;

% حلقه تکرار IDPSO
for iter = 1:max_iter
    w = w_max - ((w_max - w_min) * iter / max_iter); % وزن اینرسی پویا
    
    % محاسبه تابع هدف برای هر ذره
    for i = 1:num_particles
        fitness_value = fitness_function(positions(i, :));
        if fitness_value < p_best_values(i)
            p_best_values(i) = fitness_value;
            p_best_positions(i, :) = positions(i, :);
        end
        if fitness_value < g_best_value
            g_best_value = fitness_value;
            g_best_position = positions(i, :);
        end
    end
    
    % به‌روزرسانی سرعت و موقعیت‌ها
    for i = 1:num_particles
        r1 = rand(1, 3);
        r2 = rand(1, 3);
        velocities(i, :) = w * velocities(i, :) + ...
            c1 * r1 .* (p_best_positions(i, :) - positions(i, :)) + ...
            c2 * r2 .* (g_best_position - positions(i, :));
        positions(i, :) = positions(i, :) + velocities(i, :);
        
        % محدود کردن موقعیت به محدوده جستجو
        positions(i, :) = max(positions(i, :), bounds(:,1)');
        positions(i, :) = min(positions(i, :), bounds(:,2)');
    end
    
    % نمایش نتایج میان‌دوره‌ای
    fprintf('Iteration %d: Best Fitness = %.6f\n', iter, g_best_value);
end

% نمایش نتایج نهایی
fprintf('Best Rs: %.6f Ohm\n', g_best_position(1));
fprintf('Best Rp: %.6f Ohm\n', g_best_position(2));
fprintf('Best Io: %.6e A\n', g_best_position(3));

% تابع محاسبه خطا بر اساس نقاط کلیدی
function error = calculate_error(params, Voc, Isc, VMPP, IMPP, Nc)
    Rs = params(1);
    Rp = params(2);
    Io = params(3);
    q = 1.602176634e-19; % بار الکترون
    k = 1.380649e-23;    % ثابت بولتزمن
    T = 298.15;          % دمای سلول (کلوین)
    a = 1;               % ضریب شناسایی ایده‌آل دیود

    % معادله جریان-ولتاژ برای سه نقطه کلیدی
    Ipv = Isc;
    err_oc = Ipv - Io * (exp(q * Voc / (a * Nc * k * T)) - 1) - Voc / Rp;
    err_sc = Isc - Io * (exp(q * Rs * Isc / (a * Nc * k * T)) - 1) - Rs * Isc / Rp;
    err_mpp = IMPP - Io * (exp(q * (VMPP + Rs * IMPP) / (a * Nc * k * T)) - 1) - (VMPP + Rs * IMPP) / Rp;
    
    % محاسبه مجموع خطاها
    error = err_oc^2 + err_sc^2 + err_mpp^2;
end
