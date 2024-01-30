module Integrate

using Printf

function step!(config, model, buffer)

    #getting correct index
    t = buffer.currentIndex
    if t<=0
        t = t%buffer.size+buffer.size
    end
    
    #The actual computations

    qn = buffer.data[:, t, :, :, :]
    
    if config.integrator == :SSPRK3
        q1 = qn + model.deltaT * model.RHS(qn, config) #qn + ΔtF[qn]
        q2 = 3. / 4. * qn .+ 1. /4. * q1 .+ 1. / 4. * model.deltaT * model.RHS(q1, config) #3/4qn + 1/4q1 + 1/4ΔtF[q1]

        buffer.data[:, buffer.currentIndex + 1, :,:,:] = 1. / 3. * qn .+ 2. / 3. * q2 .+ 2. / 3. * model.deltaT * model.RHS(q2, config) #1/3qn + 2/3q2 + 2/3ΔtF[q2]
        
    elseif config.integrator == :SSPRK2

        q1 = qn .+ model.deltaT * model.RHS(qn, config) #qn +ΔtF[qn]

        buffer.data[:, buffer.currentIndex + 1, :,:,:] = 1. / 2. * qn .+ 1. / 2. * q1 .+ 1. /2. * model.deltaT * model.RHS(q1, config) # 1/2qn + 1/2q1 + 1/2ΔtF[q1]
    
    end
        
end

end