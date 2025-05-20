var net = require('net');

net.createServer(function(socket){
    socket.on('data', function(data){
        socket.write(data.toString())
    });
}).listen(55555);
