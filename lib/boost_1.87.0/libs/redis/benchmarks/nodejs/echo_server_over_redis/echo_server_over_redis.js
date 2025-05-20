import { createClient } from 'redis';
import * as net from 'net';

const client = createClient({url: 'redis://aedis.occase.de:63799' });
client.on('error', (err) => console.log('Redis Client Error', err));
await client.connect();

net.createServer(function(socket){
   socket.on('data', async function(data) {
      const value = await client.ping(data.toString());
      socket.write(data)
   });
}).listen(55555);
